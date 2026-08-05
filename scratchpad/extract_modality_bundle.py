#!/usr/bin/env python3
"""Extract the 8 sc modality feature blocks + labels + donor groups + holdout into ONE portable bundle,
so the per-modality base-learner sweep can run standalone on a desktop (no atlas, no cluster).
Mirrors train_predictor.load_block exactly."""
import os, sys, pickle, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline"))
import numpy as np, pandas as pd
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics, udon_features as UF

OUT = r"C:\Users\krog5w\.gemini\antigravity\scratch\aml-modality-bakeoff\modality_bundle.pkl"
os.makedirs(os.path.dirname(OUT), exist_ok=True)
MODS = ["RNA", "Composition", "ADT", "Lipid", "Metabolite", "GRN", "LSC", "Cell-comm"]
t0 = time.time()
ctx = build_context(Config(run_id="single_modality"))
samples = ctx.tables["samples"]; dg = samples["donor_group"].astype(str); hold = set(ctx.holdout)

def load_block(mod):
    if mod == "RNA": return UF.canonical_rna(ctx)
    c = ctx.path("_sl_%s.pkl" % mod)
    if os.path.exists(c): return pd.read_pickle(c)
    if mod == "Composition": return D._sample_level_matrix(ctx, "composition", set(samples.index))
    if mod in ("ADT", "GRN"): return dataio.sample_modality_matrix(ctx, mod)
    if mod in ("Lipid", "Metabolite"): return dataio.sample_modality_matrix(ctx, mod, min_spearman=0.3)
    if mod == "Cell-comm": return dataio.cellcomm_matrix(ctx)
    if mod == "LSC":
        t = ctx.tables.get("lsc_calls"); cols = [c for c in ["Prob_m-LSC", "Prob_p+m-LSC", "Prob_p-LSC", "MaxProb"] if c in t.columns]
        return t[cols].apply(pd.to_numeric, errors="coerce")

VARCAP = 8000                                                     # cap ultra-high-dim blocks (Cell-comm 141k) by variance
BLK = {}
for m in MODS:
    try:
        b = load_block(m).fillna(0.0); b = b[~b.index.duplicated(keep="first")].astype(np.float32)
        if b.shape[1] > VARCAP:                                    # keep top-VARCAP most-variable cols (diff_select(500) picks from these anyway)
            top = b.var(0).sort_values(ascending=False).index[:VARCAP]
            b = b[top]; print("  %-11s %s  (capped from wide block)" % (m, b.shape))
        else:
            print("  %-11s %s" % (m, b.shape))
        BLK[m] = b
    except Exception as e:
        print("  skip %s: %s" % (m, e))
MODS = [m for m in MODS if m in BLK]

M = ctx.tables.get("mutations") or genetics.build_mutation_matrix(ctx)
_m01 = {"present": 1.0, "absent": 0.0}
allidx = set(BLK["RNA"].index)
LAB = {}
for f in sorted(c for c in M.columns if str(c).startswith(("mut_", "cyto_"))):
    yr = D._labels_for_field_raw(ctx, f).map(_m01)                 # all samples (incl holdout), 1/0/NaN
    tr = [s for s in allidx if pd.notna(yr.get(s)) and s not in hold]
    yv = np.array([int(yr[s]) for s in tr])
    if (yv == 1).sum() >= 8 and (yv == 0).sum() >= 8:              # deployed well-powered floor
        LAB[f.replace("mut_", "").replace("cyto_", "")] = yr
labels = pd.DataFrame(LAB)
print("trainable mutations: %d" % labels.shape[1])

bundle = {"mods": MODS, "blocks": {m: BLK[m] for m in MODS},
          "labels": labels, "donor_group": dg, "holdout": sorted(hold),
          "all_index": sorted(allidx)}
with open(OUT, "wb") as fh:
    pickle.dump(bundle, fh, protocol=4)
print("wrote %s  (%.1f MB, %.0fs)" % (OUT, os.path.getsize(OUT) / 1e6, time.time() - t0))
