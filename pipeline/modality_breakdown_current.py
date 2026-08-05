#!/usr/bin/env python3
"""Standalone per-(driver x modality) AUROC on the CURRENT model and labels.

The earlier breakdown predates the NYU-2 label correction and the BeatAML augmentation, so its numbers
no longer match the deployed model. This recomputes, for every driver and every modality block, the
donor-grouped CV-OOF AUROC of that modality ALONE — the quantity that tells you which measurement type
actually carries a given lesion — alongside the fusion weight the deployed model assigns it.

Recipe is identical to the deployed per-modality path: StandardScaler -> differential top-500 ->
LinearSVC(C=0.02, balanced) -> percentile -> donor-grouped 3-fold CV.

  bsub -q test -W 400 -M 64000 -R "rusage[mem=64000]" -o mb.log \
    /usr/local/anaconda3-2020/bin/python modality_breakdown_current.py
"""
import os, sys, json, warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics, udon_features as UF
from amlmm.predictor import diff_select, _pct
def log(m): print(m, flush=True)

OUT = os.path.join(ROOT, "deliverables", "modality_breakdown_current.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
MODS = ["RNA", "Composition", "ADT", "Lipid", "Metabolite", "GRN", "LSC", "Cell-comm", "BulkRNA"]

ctx = build_context(Config(run_id="single_modality"))
samples = ctx.tables["samples"]; dg = samples["donor_group"].astype(str); hold = set(ctx.holdout)
import bulk_features as BF

def load_block(mod):
    if mod == "RNA": return UF.canonical_rna(ctx)
    c = ctx.path("_sl_%s.pkl" % mod)
    if os.path.exists(c): return pd.read_pickle(c)
    if mod == "BulkRNA":
        b = BF.bulk_rna_matrix(ctx); b.to_pickle(c); return b
    if mod == "Composition": return D._sample_level_matrix(ctx, "composition", set(samples.index))
    if mod in ("ADT", "GRN"): return dataio.sample_modality_matrix(ctx, mod)
    if mod in ("Lipid", "Metabolite"): return dataio.sample_modality_matrix(ctx, mod, min_spearman=0.3)
    if mod == "Cell-comm": return dataio.cellcomm_matrix(ctx)
    if mod == "LSC":
        t = ctx.tables.get("lsc_calls")
        cols = [c for c in ["Prob_m-LSC","Prob_p+m-LSC","Prob_p-LSC","MaxProb"] if c in t.columns]
        return t[cols].apply(pd.to_numeric, errors="coerce")

BLK = {}
for m in MODS:
    try:
        b = load_block(m).fillna(0.0)
        BLK[m] = b[~b.index.duplicated(keep="first")]
        log("loaded %-12s %s" % (m, BLK[m].shape))
    except Exception as e:
        log("skip %-12s %s" % (m, str(e)[:70]))
MODS = [m for m in MODS if m in BLK]

M = ctx.tables.get("mutations") or genetics.build_mutation_matrix(ctx)
_m01 = {"present": 1.0, "absent": 0.0}
def auc(y, s):
    s = np.asarray(s, float); ok = ~np.isnan(s)
    return roc_auc_score(y[ok], s[ok]) if (ok.sum() >= 4 and len(set(y[ok])) == 2) else np.nan

def cv_auc(B, ids, y, grp):
    X = B.loc[ids].values; keep = X.std(0) > 0
    if keep.sum() < 2 or len(set(grp)) < 2: return None
    oof = np.full(len(ids), np.nan)
    for tri, vai in GroupKFold(min(3, len(set(grp)))).split(X, y, grp):
        if len(set(y[tri])) < 2: continue
        sc = StandardScaler().fit(X[tri][:, keep])
        Ztr = sc.transform(X[tri][:, keep]); Zva = sc.transform(X[vai][:, keep])
        sel = diff_select(Ztr, y[tri], 500)
        oof[vai] = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(
            Ztr[:, sel], y[tri]).decision_function(Zva[:, sel])
    ok = ~np.isnan(oof)
    if ok.sum() < 4: return None
    p = _pct(np.sort(oof[ok]), oof[ok])
    a = auc(y[ok], p)
    return float(max(a, 1 - a)) if a == a else None

MUTS = []
for f in sorted(c for c in M.columns if str(c).startswith(("mut_", "cyto_"))):
    ym = D.labels_for_field(ctx, f).map(_m01)
    tr = [s for s in BLK["RNA"].index if pd.notna(ym.get(s)) and s not in hold]
    yv = np.array([int(ym[s]) for s in tr])
    if (yv == 1).sum() >= 8 and (yv == 0).sum() >= 8: MUTS.append(f)
log("drivers: %d | modalities: %s" % (len(MUTS), MODS))

# fusion weights from the deployed model, for cross-reference
try:
    PW = {m: (a.get("fused_all") or {}).get("weights", {})
          for m, a in json.load(open(os.path.join(ROOT, "deliverables",
                                "production_fused_model.json")))["per_mutation"].items()}
except Exception:
    PW = {}

res = {}
for mflag in MUTS:
    short = mflag.replace("mut_", "").replace("cyto_", "")
    yall = D._labels_for_field_raw(ctx, mflag).map(_m01)
    ym = D.labels_for_field(ctx, mflag).map(_m01)
    row = {}
    for mod in MODS:
        B = BLK[mod]
        ids = [s for s in B.index if pd.notna(ym.get(s)) and s not in hold]
        y = np.array([int(yall[s]) for s in ids])
        if (y == 1).sum() < 5 or (y == 0).sum() < 5: continue
        a = cv_auc(B, ids, y, dg.loc[ids].values)
        if a is not None: row[mod] = round(a, 4)
    if not row: continue
    best = max(row, key=row.get)
    npos = int(sum(1 for s in BLK["RNA"].index if s not in hold and yall.get(s) == 1))
    res[short] = {"standalone_auroc": row, "best_modality": best, "best_single_auroc": row[best],
                  "n_pos": npos, "fusion_weights": PW.get(short, {}),
                  "spread": round(max(row.values()) - min(row.values()), 4)}
    log("  %-14s best=%-12s %.3f | spread %.3f | %s" %
        (short, best, row[best], res[short]["spread"],
         " ".join("%s:%.2f" % (k, v) for k, v in sorted(row.items(), key=lambda kv: -kv[1])[:4])))

json.dump({"generated": time.strftime("%Y-%m-%d %H:%M"), "modalities": MODS,
           "note": "Standalone donor-grouped CV-OOF AUROC per (driver x modality) on the current "
                   "labels (NYU-2 corrected). fusion_weights are from the deployed augmented model.",
           "drivers": res}, open(OUT, "w"), indent=1)
log("\nwrote %s (%d drivers)" % (OUT, len(res)))
log("MODALITY BREAKDOWN OK")
