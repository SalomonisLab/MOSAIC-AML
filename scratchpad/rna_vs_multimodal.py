#!/usr/bin/env python3
"""Is the multimodal machinery earning its keep over RNA alone?

Three trained systems, SAME 26 mutations, SAME recipe (per-modality StandardScaler -> diff_select(500)
-> LinearSVC(C=0.02) -> percentile -> donor-grouped CV-OOF -> ridge-NNLS fusion), differing only in which
modality blocks they were given:
    model_card_bulk.json     modalities = [BulkRNA]                      <- bulk RNA-seq ALONE
    model_card.json          modalities = the 8 sc-derived blocks        <- WHAT IS DEPLOYED
    model_card_bulk_sc.json  modalities = BulkRNA + all 8                <- everything
So a paired per-mutation comparison is exact, not an apples-to-oranges guess.

Also pulls the deployed system's OWN internals (deployed_modality_breakdown.csv) to ask the same
question a second way: the sc RNA block alone vs the fused combiner, within one model.

train_auc  = donor-grouped CV-OOF AUROC on the training cohort  (the reliable metric)
heldout_auc= the same combiner on the sealed held-out set       (tiny: 26-29 samples, few positives -> noisy)
"""
import os, json
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPE = os.path.join(ROOT, "pipeline")


def card(name):
    c = json.load(open(os.path.join(PIPE, name)))
    return {m: v for m, v in c["mutations"].items()}, (c.get("meta") or {}).get("modalities")


bulk, mb = card("model_card_bulk.json")
allm, ma = card("model_card.json")
both, mo = card("model_card_bulk_sc.json")
print("BulkRNA alone   :", mb)
print("DEPLOYED (all)  :", ma or "(8 sc blocks: RNA, Composition, ADT, Lipid, Metabolite, GRN, LSC, Cell-comm)")
print("BulkRNA + all   :", mo)
print()

muts = sorted(set(bulk) & set(allm))
rows = []
for m in muts:
    rows.append(dict(
        mutation=m, n_pos=allm[m].get("n_pos"),
        bulk_tr=bulk[m].get("train_auc"), all_tr=allm[m].get("train_auc"), both_tr=both.get(m, {}).get("train_auc"),
        bulk_ho=bulk[m].get("heldout_auc"), all_ho=allm[m].get("heldout_auc"),
    ))
df = pd.DataFrame(rows)

print("=" * 92)
print("%-22s %5s | %8s %8s %8s | %8s %8s" %
      ("mutation", "n+", "BULK-only", "ALL(dep)", "diff", "BULKho", "ALLho"))
print("=" * 92)
d = df.dropna(subset=["bulk_tr", "all_tr"]).copy()
d["diff"] = d["all_tr"] - d["bulk_tr"]
for _, r in d.sort_values("diff").iterrows():
    print("%-22s %5s | %8.3f %8.3f %+8.3f | %8s %8s" %
          (r["mutation"], r["n_pos"], r["bulk_tr"], r["all_tr"], r["diff"],
           ("%.3f" % r["bulk_ho"]) if pd.notna(r["bulk_ho"]) else "-",
           ("%.3f" % r["all_ho"]) if pd.notna(r["all_ho"]) else "-"))

print()
print("=" * 60)
print("TRAIN CV-OOF (the reliable metric), paired over %d mutations" % len(d))
print("=" * 60)
print("  BulkRNA alone      : %.3f" % d["bulk_tr"].mean())
print("  ALL 8 (deployed)   : %.3f" % d["all_tr"].mean())
print("  difference         : %+.3f  (deployed - bulk-only)" % (d["all_tr"].mean() - d["bulk_tr"].mean()))
print("  deployed wins      : %d/%d   bulk-only wins: %d" %
      (int((d["diff"] > 0).sum()), len(d), int((d["diff"] < 0).sum())))
try:
    from scipy.stats import wilcoxon
    st, p = wilcoxon(d["all_tr"], d["bulk_tr"])
    print("  Wilcoxon signed-rank p = %.4g  %s" % (p, "(significant)" if p < 0.05 else "(NOT significant)"))
except Exception as e:
    print("  (wilcoxon unavailable: %s)" % e)

b = df.dropna(subset=["both_tr", "all_tr"])
if len(b):
    print()
    print("  adding BulkRNA on TOP of the 8: %.3f vs %.3f deployed  (%+.3f)"
          % (b["both_tr"].mean(), b["all_tr"].mean(), b["both_tr"].mean() - b["all_tr"].mean()))

h = df.dropna(subset=["bulk_ho", "all_ho"])
if len(h):
    print()
    print("HELD-OUT (sealed, but tiny -> noisy), paired over %d mutations" % len(h))
    print("  BulkRNA alone: %.3f | ALL 8 (deployed): %.3f  (%+.3f)"
          % (h["bulk_ho"].mean(), h["all_ho"].mean(), h["all_ho"].mean() - h["bulk_ho"].mean()))
    print("  deployed wins %d/%d" % (int((h["all_ho"] > h["bulk_ho"]).sum()), len(h)))

# ---- second angle: inside the deployed model, sc RNA block alone vs the fused combiner ----
p = os.path.join(ROOT, "scratchpad", "deployed_modality_breakdown.csv")
if os.path.exists(p):
    bd = pd.read_csv(p)
    q = bd.dropna(subset=["RNA", "fused_train_auc"])
    print()
    print("=" * 60)
    print("SECOND ANGLE: inside the deployed model (same CV, %d mutations)" % len(q))
    print("=" * 60)
    print("  sc RNA block ALONE : %.3f" % q["RNA"].mean())
    print("  fused (all 8)      : %.3f" % q["fused_train_auc"].mean())
    print("  difference         : %+.3f" % (q["fused_train_auc"].mean() - q["RNA"].mean()))
    print("  fused beats RNA-alone in %d/%d mutations" %
          (int((q["fused_train_auc"] > q["RNA"]).sum()), len(q)))
    print("  best SINGLE modality mean: %.3f  -> fusion adds %+.3f over the best single"
          % (q["best_single_auc"].mean(), q["fused_train_auc"].mean() - q["best_single_auc"].mean()))
    print()
    print("  per-modality mean CV-OOF AUROC (how good is each block on its own):")
    for mod in ["RNA", "Composition", "ADT", "Lipid", "Metabolite", "GRN", "LSC", "Cell-comm"]:
        if mod in bd.columns:
            print("     %-11s %.3f   (best modality for %d/%d mutations)"
                  % (mod, bd[mod].mean(), int((bd["best_modality"] == mod).sum()), len(bd)))
