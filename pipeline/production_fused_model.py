#!/usr/bin/env python3
"""PRODUCTION fused model: all atlas modalities + BeatAML-augmented shared blocks, augmenting SELECTIVELY.

Motivation (from the three preceding experiments):
  * BeatAML alone -> single-cell transfers poorly (AUROC 0.66) — it cannot replace single-cell training.
  * Adding BeatAML to the shared blocks helps a lot where single-cell positives are SCARCE
    (WT1 +0.226, CEBPA +0.095, IDH1 +0.063 ...) and slightly HURTS mutations the atlas already fits
    well (TP53 -0.041, TET2 -0.070) — bulk dilutes a well-fit single-cell signal.
  => augment per-mutation, only when the atlas is thin. The rule uses TRAINING label counts only
     (n_pos_atlas < AUG_MAX_POS), never test performance, so it is not selection bias.

ARMS (identical samples + donor-grouped folds):
  deployed          the 8 atlas modalities, atlas-trained          <- current production baseline
  fused_all         + BulkRNA block, BeatAML-augmented on EVERY shared block
  fused_selective   + BulkRNA block, BeatAML-augmented ONLY when n_pos_atlas < AUG_MAX_POS

Shared (poolable) blocks: BulkRNA (gene symbols), ADT (the rna2adt "Hu." prefix is stripped to match
the atlas naming), Metabolite, GRN. Lipid cannot be pooled: the atlas was imputed with the
`newnormelastic` bundle which is absent here, so its species names do not correspond to `otherelastic`.

  bsub -q test -W 700 -M 72000 -R "rusage[mem=72000]" -o prod.log \
    /usr/local/anaconda3-2020/bin/python production_fused_model.py
"""
import os, sys, json, warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
COMP = os.path.join(ROOT, "engine-code", "altanalyze3", "altanalyze3", "components")
ALT = os.path.join(ROOT, "engine-code", "altanalyze3")
sys.path.insert(0, ALT); sys.path.insert(0, COMP); sys.path.insert(0, HERE)
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics, udon_features as UF
from amlmm.predictor import diff_select, _pct
from amlmm.bulk_predictor import BulkMutationPredictor
def log(m): print(m, flush=True)

OUT = os.path.join(ROOT, "deliverables", "production_fused_model.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
AUG_MAX_POS = int(os.environ.get("AMLMM_AUG_MAX_POS", "25"))     # augment only mutations thinner than this
ATLAS8 = ["RNA", "Composition", "ADT", "Lipid", "Metabolite", "GRN", "LSC", "Cell-comm"]
POOLABLE = ["BulkRNA", "ADT", "Metabolite", "GRN"]

ctx = build_context(Config(run_id="single_modality"))
samples = ctx.tables["samples"]; dg = samples["donor_group"].astype(str); hold = set(ctx.holdout)
import bulk_features as BF

def atlas_block(mod):
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
        cols = [c for c in ["Prob_m-LSC", "Prob_p+m-LSC", "Prob_p-LSC", "MaxProb"] if c in t.columns]
        return t[cols].apply(pd.to_numeric, errors="coerce")

ATL = {}
for m in ATLAS8 + ["BulkRNA"]:
    try:
        b = atlas_block(m).fillna(0.0)
        b = b[~b.index.duplicated(keep="first")]; b = b.loc[:, ~b.columns.duplicated()]
        ATL[m] = b; log("atlas %-11s %s" % (m, b.shape))
    except Exception as e:
        log("atlas %-11s FAILED %s" % (m, str(e)[:90]))

# ---------------- BeatAML side ----------------
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
e2s = {v: k for k, v in (BP.sym2ens or {}).items()}
ba_s = [str(s) for s in d["ba_samples"]]; gs = [str(g) for g in d["genes"]]
sy = [e2s.get(g) for g in gs]; kp = [i for i, s in enumerate(sy) if s]
BA_SYM = pd.DataFrame(d["ba_X"][:, kp].astype(float), index=ba_s, columns=[sy[i] for i in kp])
BA_SYM = BA_SYM.T.groupby(level=0).sum().T
from altanalyze3.components.rna2adt.api import Rna2AdtBundle
from altanalyze3.components.rna2grn.api import Rna2GrnBundle
from altanalyze3.components.rna2metabolite.api import load_bundle as load_metab
BAB = {"BulkRNA": BA_SYM}
for name, loader in [("ADT", lambda: Rna2AdtBundle.load(os.path.join(COMP, "rna2adt", "rna2adt_bm_bundle.pkl"))),
                     ("GRN", lambda: Rna2GrnBundle.load(os.path.join(COMP, "rna2grn", "rna2grn_bundle.pkl.gz"))),
                     ("Metabolite", lambda: load_metab(os.path.join(COMP, "rna2metabolite", "artifacts", "rna2metabolite_aml_bundle.pkl.gz")))]:
    try:
        r = loader().predict_from_dataframe(BA_SYM)
        df = (r.predictions if hasattr(r, "predictions") else r).reindex(BA_SYM.index).fillna(0.0)
        BAB[name] = df.loc[:, ~df.columns.duplicated()]; log("beataml %-11s %s" % (name, BAB[name].shape))
    except Exception as e:
        log("beataml %-11s FAILED %s" % (name, str(e)[:90]))

def nkey(s):
    s = str(s)
    if s.startswith("Hu."): s = s[3:]          # rna2adt prefixes every antibody with "Hu."
    return "".join(ch for ch in s.upper() if ch.isalnum())
SHARED, BCOL = {}, {}
for m in POOLABLE:
    if m not in ATL or m not in BAB: continue
    bmap = {}
    for c in BAB[m].columns: bmap.setdefault(nkey(c), c)
    pairs, seen = [], set()
    for c in ATL[m].columns:
        k = nkey(c)
        if k in bmap and k not in seen:
            seen.add(k); pairs.append((c, bmap[k]))
    if len(pairs) >= 20:
        SHARED[m] = [a for a, _ in pairs]; BCOL[m] = [b for _, b in pairs]
        log("POOLABLE %-11s %d shared features" % (m, len(pairs)))
    else:
        log("POOLABLE %-11s only %d shared -> not pooled" % (m, len(pairs)))

# ---------------- labels ----------------
M = ctx.tables.get("mutations") or genetics.build_mutation_matrix(ctx)
_m01 = {"present": 1.0, "absent": 0.0}
cats = [str(c) for c in d["drivers"]]; baL = pd.DataFrame(d["ba_L"].astype(float), index=ba_s, columns=cats)
def gene_of(c):
    cl = str(c).lower()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl: return "INV16"
    if "kmt2a" in cl: return "KMT2A"
    return str(c).split("_")[0].split("-")[0].upper()
def ba_labels_for(short):
    s = short.upper().replace("-", "").replace("_", "")
    cc = [c for c in cats if gene_of(c).replace("-", "").replace("_", "") == s]
    if short.upper().startswith("FLT3-ITD"): cc = [c for c in cats if "ITD" in c.upper()]
    if short.upper().startswith("FLT3-TKD"): cc = [c for c in cats if "TKD" in c.upper()]
    if not cc: return None
    sub = baL[cc]; y = (sub.max(axis=1) == 1).astype(float); y[sub.isna().all(axis=1)] = np.nan
    return y
def auc(y, s):
    s = np.asarray(s, float); ok = ~np.isnan(s)
    return roc_auc_score(y[ok], s[ok]) if (ok.sum() >= 4 and len(set(y[ok])) == 2) else np.nan
def conf(y, call):
    y = np.asarray(y).astype(int); call = np.asarray(call).astype(int)
    tp = int(((call==1)&(y==1)).sum()); fp = int(((call==1)&(y==0)).sum())
    fn = int(((call==0)&(y==1)).sum()); tn = int(((call==0)&(y==0)).sum())
    se = tp/(tp+fn) if tp+fn else 0.0; sp = tn/(tn+fp) if tn+fp else 0.0
    pr = tp/(tp+fp) if tp+fp else 0.0
    return dict(tp=tp,fp=fp,fn=fn,tn=tn,sensitivity=round(se,4),specificity=round(sp,4),
                precision=round(pr,4), f1=round(2*pr*se/(pr+se),4) if pr+se else 0.0)
def nested_thr_metrics(comb, y, grp):
    ng = min(3, len(set(grp)))
    if ng < 2: return conf(y, comb >= 0.5)
    call = np.zeros(len(y), int)
    for tri, tei in GroupKFold(ng).split(comb, y, grp):
        bt, bf = 0.5, -1.0
        for t in np.unique(comb[tri]):
            c = conf(y[tri], comb[tri] >= t)
            if c["f1"] > bf: bf, bt = c["f1"], float(t)
        call[tei] = (comb[tei] >= bt).astype(int)
    return conf(y, call)

MUTS = []
for f in sorted(c for c in M.columns if str(c).startswith(("mut_", "cyto_"))):
    ym = D.labels_for_field(ctx, f).map(_m01)
    tr = [s for s in ATL["RNA"].index if pd.notna(ym.get(s)) and s not in hold]
    yv = np.array([int(ym[s]) for s in tr])
    if (yv == 1).sum() >= 8 and (yv == 0).sum() >= 8: MUTS.append(f)
log("trainable mutations: %d | poolable: %s | AUG_MAX_POS=%d" % (len(MUTS), list(SHARED), AUG_MAX_POS))

ARMS = ["deployed", "fused_all", "fused_selective"]
res = {a: {} for a in ARMS}; augmented_for = {}
for mflag in MUTS:
    short = mflag.replace("mut_", "").replace("cyto_", "")
    yall = D._labels_for_field_raw(ctx, mflag).map(_m01); ym = D.labels_for_field(ctx, mflag).map(_m01)
    yba = ba_labels_for(short); nba = int(np.nansum(yba.values)) if yba is not None else 0
    npos_atlas = int(sum(1 for s in ATL["RNA"].index if s not in hold and yall.get(s) == 1))
    for arm in ARMS:
        blocks = ATLAS8 if arm == "deployed" else (ATLAS8 + ["BulkRNA"])
        do_aug = (arm != "deployed") and yba is not None and nba >= 5
        if arm == "fused_selective":
            do_aug = do_aug and (npos_atlas < AUG_MAX_POS)
        if arm == "fused_selective": augmented_for[short] = bool(do_aug)
        oof_p = {}
        for mod in blocks:
            if mod not in ATL: continue
            poolable = do_aug and (mod in SHARED)
            cols = SHARED[mod] if poolable else list(ATL[mod].columns)
            A = ATL[mod][cols] if poolable else ATL[mod]
            ids = [s for s in A.index if pd.notna(ym.get(s)) and s not in hold]
            y = np.array([int(yall[s]) for s in ids])
            if (y == 1).sum() < 5 or (y == 0).sum() < 5: continue
            grp = dg.loc[ids].values; Xa = A.loc[ids].values
            if poolable:
                bsel = yba.dropna().index
                Xb = BAB[mod].loc[bsel, BCOL[mod]].values; yb = yba.loc[bsel].values.astype(int)
                mb, sb = Xb.mean(0), Xb.std(0); sb[sb == 0] = 1.0
                Zb = (Xb - mb) / sb
            oof = np.full(len(ids), np.nan)
            ng = min(3, len(set(grp)))
            if ng < 2: continue
            for tri, vai in GroupKFold(ng).split(Xa, y, grp):
                if len(set(y[tri])) < 2: continue
                sc = StandardScaler().fit(Xa[tri]); Ztr = sc.transform(Xa[tri]); Zva = sc.transform(Xa[vai])
                Xtr, ytr = Ztr, y[tri]
                if poolable:
                    Xtr = np.vstack([Ztr, Zb]); ytr = np.concatenate([y[tri], yb])
                k2 = Xtr.std(0) > 0
                if k2.sum() < 2: continue
                sel = diff_select(Xtr[:, k2], ytr, 500)
                est = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(Xtr[:, k2][:, sel], ytr)
                oof[vai] = est.decision_function(Zva[:, k2][:, sel])
            ok = ~np.isnan(oof)
            if ok.sum() < 4: continue
            p = np.full(len(ids), np.nan); p[ok] = _pct(np.sort(oof[ok]), oof[ok])
            a0 = auc(y[ok], p[ok])
            if a0 == a0 and a0 < 0.5: p = 1 - p
            oof_p[mod] = dict(zip(ids, p))
        if not oof_p: continue
        common = sorted(set.intersection(*[set(oof_p[m]) for m in oof_p]))
        yco = np.array([int(yall[s]) for s in common])
        if len(common) < 8 or len(set(yco)) != 2: continue
        gco = dg.loc[common].values
        cols_ = list(oof_p); O = np.column_stack([[oof_p[m][s] for s in common] for m in cols_])
        Am = np.vstack([O, np.eye(O.shape[1])]); b = np.concatenate([yco.astype(float), np.zeros(O.shape[1])])
        w, _ = nnls(Am, b); w = w / (w.sum() or 1.0); comb = O @ w
        a1 = auc(yco, comb); nm = nested_thr_metrics(comb, yco, gco)
        nm.update(mutation=short, auroc=round(float(a1), 4) if a1 == a1 else None,
                  n_pos_atlas=npos_atlas, n_pos_beataml=nba, augmented=bool(do_aug),
                  weights={m: round(float(x), 3) for m, x in zip(cols_, w) if x > 0.01})
        res[arm][short] = nm
    log("  %-14s n+atlas=%-3d n+beat=%-4d aug=%-5s | %s" % (short, npos_atlas, nba,
        augmented_for.get(short), " ".join("%s=%.3f" % (a[:4], res[a][short]["auroc"])
        for a in ARMS if short in res[a] and res[a][short]["auroc"] is not None)))

sh = sorted(set.intersection(*[set(res[a]) for a in ARMS]))
summ = {}
for a in ARMS:
    v = [res[a][m] for m in sh]
    summ[a] = {k: round(float(np.mean([x[k] for x in v if x.get(k) is not None])), 4)
               for k in ("auroc", "sensitivity", "specificity", "f1")}
    tp = sum(x["tp"] for x in v); fp = sum(x["fp"] for x in v); fn = sum(x["fn"] for x in v); tn = sum(x["tn"] for x in v)
    summ[a].update(n_mutations=len(sh), pooled_sensitivity=round(tp/(tp+fn), 4) if tp+fn else None,
                   pooled_specificity=round(tn/(tn+fp), 4) if tn+fp else None)
try:
    from scipy.stats import wilcoxon
    summ["wilcoxon_selective_vs_deployed_auroc_p"] = float(wilcoxon(
        [res["fused_selective"][m]["auroc"] for m in sh], [res["deployed"][m]["auroc"] for m in sh])[1])
except Exception: pass
json.dump({"generated": time.strftime("%Y-%m-%d %H:%M"), "aug_max_pos": AUG_MAX_POS,
           "poolable_blocks": {m: len(SHARED[m]) for m in SHARED},
           "n_augmented": int(sum(1 for m in sh if augmented_for.get(m))),
           "augmented_mutations": sorted(m for m in sh if augmented_for.get(m)),
           "summary": summ, "per_mutation": {m: {a: res[a][m] for a in ARMS} for m in sh}},
          open(OUT, "w"), indent=1)
log("\n=== PRODUCTION FUSED MODEL (%d mutations, nested-CV operating point) ===" % len(sh))
for a in ARMS:
    s = summ[a]; log("  %-17s AUROC %.3f  sens %.3f  spec %.3f  F1 %.3f | pooled sens %.3f spec %.3f"
                     % (a, s["auroc"], s["sensitivity"], s["specificity"], s["f1"],
                        s["pooled_sensitivity"] or 0, s["pooled_specificity"] or 0))
log("augmented %d/%d mutations (n_pos_atlas < %d): %s" % (
    sum(1 for m in sh if augmented_for.get(m)), len(sh), AUG_MAX_POS,
    sorted(m for m in sh if augmented_for.get(m))))
log("wrote %s" % OUT)
log("PRODUCTION FUSED OK")
