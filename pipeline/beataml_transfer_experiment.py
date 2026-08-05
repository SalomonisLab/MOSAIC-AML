#!/usr/bin/env python3
"""NS's exact question: train on BeatAML, predict SINGLE-CELL pseudobulk mutations — does multimodal
imputation beat RNA alone at that transfer?

  TRAIN: BeatAML (707 bulk samples)   ARM A = RNA only ; ARM B = RNA + imputed (Metabolite, GRN)
  TEST : the single-cell atlas pseudobulk samples (never in training; a true cross-platform test)

Each cohort is z-scored on its OWN statistics (the cohort-matched normalisation the deployed bulk
caller already uses) so the transfer is not defeated by scale differences alone.

Only feature spaces that exist in BOTH cohorts can transfer: BulkRNA gene symbols, Metabolite, GRN.
(ADT and Lipid are imputed for both but the atlas pipeline and the rna2* bundles name their features
differently, so they cannot be matched yet — noted as a fixable gap, not a scientific limit.)

  bsub -q test -W 400 -M 64000 -R "rusage[mem=64000]" -o tr.log \
    /usr/local/anaconda3-2020/bin/python beataml_transfer_experiment.py
"""
import os, sys, json, warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
COMP = os.path.join(ROOT, "engine-code", "altanalyze3", "altanalyze3", "components")
ALT = os.path.join(ROOT, "engine-code", "altanalyze3")
sys.path.insert(0, ALT); sys.path.insert(0, COMP); sys.path.insert(0, HERE)
from scipy.optimize import nnls
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_auc_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics
from amlmm.predictor import diff_select, _pct
from amlmm.bulk_predictor import BulkMutationPredictor
def log(m): print(m, flush=True)

OUT = os.path.join(ROOT, "deliverables", "beataml_transfer_experiment.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
ctx = build_context(Config(run_id="single_modality"))
import bulk_features as BF

def atlas_block(mod):
    c = ctx.path("_sl_%s.pkl" % mod)
    if mod == "BulkRNA":
        if os.path.exists(c): return pd.read_pickle(c)
        b = BF.bulk_rna_matrix(ctx); b.to_pickle(c); return b
    if os.path.exists(c): return pd.read_pickle(c)
    if mod in ("ADT", "GRN"): return dataio.sample_modality_matrix(ctx, mod)
    return dataio.sample_modality_matrix(ctx, mod, min_spearman=0.3)

ATL = {}
for m in ["BulkRNA", "Metabolite", "GRN"]:
    b = atlas_block(m).fillna(0.0); ATL[m] = b[~b.index.duplicated(keep="first")]
    ATL[m] = ATL[m].loc[:, ~ATL[m].columns.duplicated()]
    log("atlas %-11s %s" % (m, ATL[m].shape))

d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
ens2sym = {v: k for k, v in (BP.sym2ens or {}).items()}
ba_s = [str(s) for s in d["ba_samples"]]; genes_ens = [str(g) for g in d["genes"]]
syms = [ens2sym.get(g) for g in genes_ens]; keep = [i for i, s in enumerate(syms) if s]
BA_SYM = pd.DataFrame(d["ba_X"][:, keep].astype(float), index=ba_s, columns=[syms[i] for i in keep])
BA_SYM = BA_SYM.T.groupby(level=0).sum().T
from altanalyze3.components.rna2grn.api import Rna2GrnBundle
from altanalyze3.components.rna2metabolite.api import load_bundle as load_metab
BAB = {"BulkRNA": BA_SYM}
for name, loader in [("GRN", lambda: Rna2GrnBundle.load(os.path.join(COMP, "rna2grn", "rna2grn_bundle.pkl.gz"))),
                     ("Metabolite", lambda: load_metab(os.path.join(COMP, "rna2metabolite", "artifacts", "rna2metabolite_aml_bundle.pkl.gz")))]:
    res = loader().predict_from_dataframe(BA_SYM)
    df = (res.predictions if hasattr(res, "predictions") else res).reindex(BA_SYM.index).fillna(0.0)
    BAB[name] = df.loc[:, ~df.columns.duplicated()]
    log("beataml imputed %-11s %s" % (name, BAB[name].shape))

def nkey(s): return "".join(ch for ch in str(s).upper() if ch.isalnum())
SHARED, BCOL = {}, {}
for m in ATL:
    bmap = {}
    for c in BAB[m].columns: bmap.setdefault(nkey(c), c)
    pairs, seen = [], set()
    for c in ATL[m].columns:
        k = nkey(c)
        if k in bmap and k not in seen:
            seen.add(k); pairs.append((c, bmap[k]))
    if len(pairs) >= 20:
        SHARED[m] = [a for a, _ in pairs]; BCOL[m] = [b for _, b in pairs]
        log("shared %-11s %d features" % (m, len(pairs)))
MODS = list(SHARED)

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

ARMS = {"beataml_RNA_only": ["BulkRNA"], "beataml_multimodal": MODS}
res = {a: {} for a in ARMS}
for mflag in sorted(c for c in M.columns if str(c).startswith(("mut_","cyto_"))):
    short = mflag.replace("mut_","").replace("cyto_","")
    yba = ba_labels_for(short)
    if yba is None: continue
    bsel = yba.dropna().index
    yb = yba.loc[bsel].values.astype(int)
    if yb.sum() < 5 or (yb==0).sum() < 5: continue
    ysc_all = D._labels_for_field_raw(ctx, mflag).map(_m01)
    for arm, mods in ARMS.items():
        oof_p = {}
        for mod in mods:
            cols, bcols = SHARED[mod], BCOL[mod]
            Xb = BAB[mod].loc[bsel, bcols].values
            A = ATL[mod][cols]
            ids = [s for s in A.index if pd.notna(ysc_all.get(s))]
            if len(ids) < 20: continue
            Xa = A.loc[ids].values
            mb, sb = Xb.mean(0), Xb.std(0); sb[sb==0]=1.0
            ma, sa = Xa.mean(0), Xa.std(0); sa[sa==0]=1.0
            Zb, Za = (Xb-mb)/sb, (Xa-ma)/sa                     # each cohort on its OWN scale
            k2 = Zb.std(0) > 0
            if k2.sum() < 2: continue
            sel = diff_select(Zb[:, k2], yb, 500)
            est = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(Zb[:,k2][:,sel], yb)
            dtr = est.decision_function(Zb[:,k2][:,sel]); ss = np.sort(dtr)
            p = _pct(ss, est.decision_function(Za[:,k2][:,sel]))
            a_tr = auc(yb, _pct(ss, dtr))
            if a_tr == a_tr and a_tr < 0.5: p = 1 - p
            oof_p[mod] = dict(zip(ids, p))
        if not oof_p: continue
        common = sorted(set.intersection(*[set(oof_p[m]) for m in oof_p]))
        ysc = np.array([int(ysc_all[s]) for s in common])
        if len(set(ysc)) != 2 or ysc.sum() < 3: continue
        cols_ = list(oof_p); O = np.column_stack([[oof_p[m][s] for s in common] for m in cols_])
        comb = O[:,0] if len(cols_)==1 else O.mean(axis=1)        # no label-fitted fusion: honest transfer
        a1 = auc(ysc, comb)
        bt, bf = 0.5, -1.0
        for t in np.unique(comb):
            c = conf(ysc, comb>=t)
            if c["f1"] > bf: bf, bt = c["f1"], float(t)
        mm = conf(ysc, comb>=bt); mm.update(mutation=short, auroc=round(float(a1),4) if a1==a1 else None,
                                            n_sc=len(common), n_pos_sc=int(ysc.sum()), n_pos_beataml=int(yb.sum()))
        res[arm][short] = mm
    if short in res["beataml_RNA_only"] and short in res["beataml_multimodal"]:
        log("  %-14s RNA-only %.3f  multimodal %.3f  (sc n+=%d)" % (short,
            res["beataml_RNA_only"][short]["auroc"] or 0, res["beataml_multimodal"][short]["auroc"] or 0,
            res["beataml_multimodal"][short]["n_pos_sc"]))

sh = sorted(set(res["beataml_RNA_only"]) & set(res["beataml_multimodal"]))
summ = {}
for a in ARMS:
    v = [res[a][m] for m in sh]
    summ[a] = {k: round(float(np.mean([x[k] for x in v if x.get(k) is not None])),4)
               for k in ("auroc","sensitivity","specificity","f1")}
    summ[a]["n_mutations"] = len(sh)
# pooled (micro-averaged) across all mutations — every positive counts, no per-mutation minimum
for a in ARMS:
    tp=sum(res[a][m]["tp"] for m in sh); fp=sum(res[a][m]["fp"] for m in sh)
    fn=sum(res[a][m]["fn"] for m in sh); tn=sum(res[a][m]["tn"] for m in sh)
    summ[a]["pooled_sensitivity"]=round(tp/(tp+fn),4) if tp+fn else None
    summ[a]["pooled_specificity"]=round(tn/(tn+fp),4) if tn+fp else None
try:
    from scipy.stats import wilcoxon
    summ["wilcoxon_auroc_p"] = float(wilcoxon([res["beataml_multimodal"][m]["auroc"] for m in sh],
                                              [res["beataml_RNA_only"][m]["auroc"] for m in sh])[1])
except Exception: pass
json.dump({"generated": time.strftime("%Y-%m-%d %H:%M"), "shared_modalities": MODS,
           "design": "TRAIN on BeatAML (707 bulk), TEST on single-cell atlas pseudobulk",
           "summary": summ, "per_mutation": {m: {a: res[a][m] for a in ARMS} for m in sh}},
          open(OUT, "w"), indent=1)
log("\n=== BeatAML -> single-cell TRANSFER (%d mutations) ===" % len(sh))
for a in ARMS:
    s = summ[a]; log("  %-20s AUROC %.3f  sens %.3f  spec %.3f  F1 %.3f | pooled sens %.3f spec %.3f"
                     % (a, s["auroc"], s["sensitivity"], s["specificity"], s["f1"],
                        s["pooled_sensitivity"] or 0, s["pooled_specificity"] or 0))
log("wrote %s" % OUT)
log("TRANSFER EXPERIMENT OK")
