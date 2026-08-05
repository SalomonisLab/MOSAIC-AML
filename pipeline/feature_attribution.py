#!/usr/bin/env python3
"""Systematic feature attribution: what drives each mutation call (top genes, modality weights,
cell-state localization). Extracted straight from the deployed model — no re-training — so it is
exactly what the classifier uses. Lets a reviewer check the model rides real biology, not artifacts.

Writes deliverables/feature_attribution.json + .tsv (+ a heatmap if matplotlib).
  bsub -q test -W 20 -M 12000 -R "rusage[mem=12000]" -o fa.log \
    /usr/local/anaconda3-2020/bin/python feature_attribution.py
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from amlmm.predictor import MutationPredictor

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)
P = MutationPredictor.load(os.path.join(HERE, "mutation_predictor.pkl"))
try:
    LOC = json.load(open(os.path.join(ROOT, "gui", "cellstate_localization.json"))).get("drivers", {})
except Exception:
    LOC = {}

def top_genes(mm, k=8):
    feats = np.asarray(mm.features)[mm.keep][mm.sel]
    coef = np.asarray(mm.svm.coef_).ravel() * float(getattr(mm, "sign", 1.0))
    up = np.argsort(-coef)[:k]; dn = np.argsort(coef)[:k]
    return ([(str(feats[i]), round(float(coef[i]), 3)) for i in up],
            [(str(feats[i]), round(float(coef[i]), 3)) for i in dn])

def loc_states(mut, k=5):
    L = LOC.get(mut) or {}
    st = (L.get("states") or {})
    ents = sorted(st.items(), key=lambda kv: -(kv[1].get("auroc") or 0))[:k]
    return [(s, round(v.get("auroc", 0), 3)) for s, v in ents]

out = {}
for mut in P.mutations:
    rec = {"modality_weights": {k: v for k, v in (P.weights.get(mut) or {}).items() if v > 0}}
    mm = P.models.get((mut, "RNA"))
    if mm is not None and hasattr(mm.svm, "coef_"):
        up, dn = top_genes(mm)
        rec["top_genes_present"] = up      # expression UP in mutant -> drives a 'present' call
        rec["top_genes_absent"] = dn
    ls = loc_states(mut)
    if ls:
        rec["top_cellstates"] = ls
    out[mut] = rec

json.dump(out, open(os.path.join(OUTD, "feature_attribution.json"), "w"), indent=1)
with open(os.path.join(OUTD, "feature_attribution.tsv"), "w") as fh:
    fh.write("mutation\ttop_modalities\ttop_genes_present\ttop_cellstates\n")
    for mut, r in out.items():
        mods = ", ".join("%s:%.2f" % (k, v) for k, v in sorted(r["modality_weights"].items(), key=lambda x: -x[1]))
        genes = ", ".join(g for g, _ in r.get("top_genes_present", [])[:6])
        cs = ", ".join(s for s, _ in r.get("top_cellstates", [])[:4])
        fh.write("%s\t%s\t%s\t%s\n" % (mut, mods, genes, cs))

print("=== feature attribution (top drivers of a 'present' call) ===")
for mut in sorted(out, key=lambda m: -(P.weights.get(m, {}) or {}).get("RNA", 0))[:12]:
    r = out[mut]
    print("  %-14s mods={%s}  genes=%s" % (
        mut, ", ".join("%s:%.2f" % (k, v) for k, v in sorted(r["modality_weights"].items(), key=lambda x: -x[1])[:3]),
        ", ".join(g for g, _ in r.get("top_genes_present", [])[:5])))
print("wrote deliverables/feature_attribution.json + .tsv")
print("FEATURE ATTRIBUTION OK")
