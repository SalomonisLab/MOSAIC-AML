#!/usr/bin/env python3
"""Nathan's experiment: impute the 4 modalities FROM BeatAML bulk RNA (rna2adt/grn/lipid/metabolite),
then does the MULTIMODAL recipe (RNA + imputed) beat RNA-alone on BeatAML?

If the imputed modalities carried information beyond the RNA they are computed from, arm B > arm A.
Prediction: ~no gain (imputed = deterministic function of the same RNA; a linear model on RNA already
sees it). Both arms use identical samples/folds/recipe, so the contrast is fair.

  bsub -q test -W 200 -M 40000 -R "rusage[mem=40000]" -o bie.log \
    /usr/local/anaconda3-2020/bin/python beataml_impute_experiment.py
"""
import os, sys, json, warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
COMP = os.path.join(ROOT, "engine-code", "altanalyze3", "altanalyze3", "components")
ALT_ROOT = os.path.join(ROOT, "engine-code", "altanalyze3")   # makes `import altanalyze3` resolvable (fixes adt/grn absolute imports)
sys.path.insert(0, ALT_ROOT); sys.path.insert(0, COMP); sys.path.insert(0, HERE)
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from amlmm.bulk_predictor import BulkMutationPredictor
from amlmm.predictor import diff_select, _pct
def log(m): print(m, flush=True)

# ---- BeatAML RNA (ENSG) -> gene symbols ----
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
ens2sym = {v: k for k, v in (BP.sym2ens or {}).items()}
genes_ens = [str(g) for g in d["genes"]]; samp = [str(s) for s in d["ba_samples"]]
syms = [ens2sym.get(g) for g in genes_ens]; keepi = [i for i, s in enumerate(syms) if s]
RNA = pd.DataFrame(d["ba_X"][:, keepi].astype(float), index=samp, columns=[syms[i] for i in keepi])
RNA = RNA.T.groupby(level=0).sum().T                          # collapse duplicate symbols
log("BeatAML RNA: %d samples x %d symbol-mapped genes (of %d ENSG); value range %.2f–%.2f mean %.2f"
    % (RNA.shape[0], RNA.shape[1], len(genes_ens), RNA.values.min(), RNA.values.max(), RNA.values.mean()))

# ---- impute the 4 modalities ----
from altanalyze3.components.rna2adt.api import Rna2AdtBundle
from altanalyze3.components.rna2grn.api import Rna2GrnBundle
from altanalyze3.components.rna2lipid.api import Rna2LipidBundle
from altanalyze3.components.rna2metabolite.api import load_bundle as load_metab

def impute():
    out = {}
    specs = [("ADT", lambda: Rna2AdtBundle.load(os.path.join(COMP, "rna2adt", "rna2adt_bm_bundle.pkl"))),
             ("GRN", lambda: Rna2GrnBundle.load(os.path.join(COMP, "rna2grn", "rna2grn_bundle.pkl.gz"))),
             ("Lipid", lambda: Rna2LipidBundle.load(os.path.join(COMP, "rna2lipid", "otherelastic_multitask_try.pkl"))),
             ("Metabolite", lambda: load_metab(os.path.join(COMP, "rna2metabolite", "artifacts", "rna2metabolite_aml_bundle.pkl.gz")))]
    for name, loader in specs:
        try:
            b = loader(); res = b.predict_from_dataframe(RNA)
            df = res.predictions if hasattr(res, "predictions") else res
            df = df.reindex(RNA.index)
            nz = float((df.abs() > 1e-9).mean().mean()); var = float(df.var().mean())
            mg = None
            try: mg = res.summary.get("matched_genes")
            except Exception: pass
            log("  imputed %-11s -> %s | matched_genes=%s | frac_nonzero=%.2f mean_var=%.4f"
                % (name, tuple(df.shape), mg, nz, var))
            if df.shape[1] >= 2 and var > 1e-8:
                out[name] = df.fillna(0.0)
            else:
                log("    DEGENERATE — dropping %s" % name)
        except Exception as e:
            log("    FAILED %s: %s" % (name, str(e)[:160]))
    return out

IMP = impute()
# CONTROL: random nonlinear RNA transforms, one per imputed modality, matched to its width. If rna_imputed
# ≈ rna_random the gain is just nonlinear feature-expansion of the SAME RNA (no biology); if rna_imputed >
# rna_random the imputers add real learned information beyond a generic nonlinear RNA transform.
rng = np.random.default_rng(0)
Xs = (RNA.values - RNA.values.mean(0)) / (RNA.values.std(0) + 1e-9)
RND = {}
for name, df in IMP.items():
    k = df.shape[1]; W = rng.standard_normal((Xs.shape[1], k)) / np.sqrt(Xs.shape[1])
    RND["rand_" + name] = pd.DataFrame(np.tanh(Xs @ W), index=RNA.index, columns=["r%d" % i for i in range(k)])
BLK = {"RNA": RNA}; BLK.update(IMP); BLK.update(RND)
imp_mods = [m for m in ("ADT", "Lipid", "Metabolite", "GRN") if m in IMP]
ARMS = {"rna_only": ["RNA"], "rna_imputed": ["RNA"] + imp_mods, "rna_random": ["RNA"] + ["rand_" + m for m in imp_mods]}
log("\nblocks: %s | arms: %s" % (list(BLK), {k: v for k, v in ARMS.items()}))

# ---- labels ----
cats = [str(c) for c in d["drivers"]]; baL = d["ba_L"].astype(float)
L = pd.DataFrame(baL, index=samp, columns=cats)

def auc(y, s):
    s = np.asarray(s, float); ok = ~np.isnan(s)
    return roc_auc_score(y[ok], s[ok]) if (ok.sum() >= 4 and len(set(y[ok])) == 2) else np.nan

def cv_oof(B, ids, y):
    X = B.loc[ids].values; keep = X.std(0) > 0
    if keep.sum() < 2: return None
    oof = np.full(len(ids), np.nan)
    for tri, vai in StratifiedKFold(5, shuffle=True, random_state=0).split(X, y):
        if len(set(y[tri])) < 2: continue
        sc = StandardScaler().fit(X[tri][:, keep]); Ztr = sc.transform(X[tri][:, keep]); Zva = sc.transform(X[vai][:, keep])
        sel = diff_select(Ztr, y[tri], 500)
        oof[vai] = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(Ztr[:, sel], y[tri]).decision_function(Zva[:, sel])
    ok = ~np.isnan(oof)
    if ok.sum() < 4: return None
    p = np.full(len(ids), np.nan); p[ok] = _pct(np.sort(oof[ok]), oof[ok])
    a = auc(y[ok], p[ok])
    if a == a and a < 0.5: p = 1 - p
    return p, (max(a, 1 - a) if a == a else 0.5)

def f1max(comb, y):
    bt, bf = 0.5, -1.0
    for t in np.unique(comb):
        pr = comb >= t; tp = int((pr & (y == 1)).sum()); fp = int((pr & (y == 0)).sum()); fn = int(((~pr) & (y == 1)).sum())
        pr_ = tp/(tp+fp) if tp+fp else 0; rc = tp/(tp+fn) if tp+fn else 0; f1 = 2*pr_*rc/(pr_+rc) if pr_+rc else 0
        if f1 > bf: bf, bt = f1, float(t)
    return bt, bf

MUTS = [c for c in cats if (L[c] == 1).sum() >= 8 and (L[c] == 0).sum() >= 8]
log("evaluable mutations: %d" % len(MUTS))
res = {}; SCORES = {}
for arm, mods in ARMS.items():
    mods = [m for m in mods if m in BLK]; per = {}
    for c in MUTS:
        y = L[c]; ids = [s for s in RNA.index if pd.notna(y.get(s))]
        yv = np.array([int(y[s]) for s in ids])
        oof_p, cvauc = {}, {}
        for m in mods:
            cv = cv_oof(BLK[m], ids, yv)
            if cv: oof_p[m] = cv[0]; cvauc[m] = cv[1]
        if not oof_p: continue
        cols = list(oof_p); O = np.column_stack([oof_p[m] for m in cols])
        if len(cols) == 1:
            comb = O[:, 0]
        else:
            A = np.vstack([O, np.eye(O.shape[1])]); b = np.concatenate([yv.astype(float), np.zeros(O.shape[1])])
            w, _ = nnls(A, b); w = w / (w.sum() or 1.0); comb = O @ w
        a = auc(yv, comb); bt, bf1 = f1max(comb, yv)
        acc = float(((comb >= bt) == (yv == 1)).mean())      # accuracy at the F1-max operating point
        per[c] = {"accuracy": round(acc, 4), "auroc": round(float(a), 4) if a == a else None,
                  "f1": round(float(bf1), 4), "n_pos": int(yv.sum()), "n": len(yv), "prevalence": round(float(yv.mean()), 4)}
        SCORES.setdefault(arm, {})[c] = {"score": [round(float(x), 5) for x in comb], "y": [int(v) for v in yv]}
    res[arm] = per

shared = sorted(set.intersection(*[set(res[a]) for a in ARMS]))
def amean(a, k):
    v = [res[a][m][k] for m in shared if res[a][m][k] is not None]
    return round(float(np.mean(v)), 4)
def dp(a1, a0, k):
    try:
        from scipy.stats import wilcoxon
        p = float(wilcoxon([res[a1][m][k] for m in shared], [res[a0][m][k] for m in shared])[1])
    except Exception:
        p = None
    return round(amean(a1, k) - amean(a0, k), 4), p
means = {a: {"accuracy": amean(a, "accuracy"), "auroc": amean(a, "auroc"), "f1": amean(a, "f1")} for a in ARMS}
d_imp = dict(zip(("delta", "p"), dp("rna_imputed", "rna_only", "accuracy")))
d_rnd = dict(zip(("delta", "p"), dp("rna_random", "rna_only", "accuracy")))
d_bio = dict(zip(("delta", "p"), dp("rna_imputed", "rna_random", "accuracy")))
out = {"generated": time.strftime("%Y-%m-%d %H:%M"), "metric": "mean per-mutation ACCURACY at the F1-max operating point",
       "imputed_modalities": list(IMP), "n_mutations": len(shared),
       "arm_means": means, "imputed_vs_only": d_imp, "random_vs_only": d_rnd,
       "imputed_vs_random_REAL_BIOLOGY": d_bio,
       "caveat": "accuracy is prevalence-dominated on rare mutations (mean prevalence ~%.1f%%), so all arms sit "
                 "near the 'call-all-absent' baseline and differences are compressed; deltas still order the arms."
                 % (100 * np.mean([res["rna_only"][m]["prevalence"] for m in shared])),
       "per_mutation": {m: {a: res[a][m] for a in ARMS} for m in shared}}
json.dump(out, open(os.path.join(ROOT, "deliverables", "beataml_impute_experiment.json"), "w"), indent=1)
json.dump(SCORES, open(os.path.join(ROOT, "deliverables", "beataml_impute_scores.json"), "w"))   # per-sample OOF -> any metric later
log("\n=== BeatAML: mean per-mutation ACCURACY over %d mutations (5-fold CV-OOF, F1-max operating point) ===" % len(shared))
for a in ARMS:
    log("  %-12s accuracy %.4f   (AUROC %.3f  F1 %.3f)" % (a, means[a]["accuracy"], means[a]["auroc"], means[a]["f1"]))
log("  imputed − only   : %+.4f  p=%s" % (d_imp["delta"], d_imp["p"]))
log("  random  − only   : %+.4f  p=%s   (nonlinear feature-expansion baseline)" % (d_rnd["delta"], d_rnd["p"]))
log("  imputed − random : %+.4f  p=%s   (== REAL biological info beyond nonlinearity)" % (d_bio["delta"], d_bio["p"]))
log("BEATAML IMPUTE EXPERIMENT OK")
