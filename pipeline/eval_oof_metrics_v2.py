#!/usr/bin/env python3
"""Provability upgrade to the single-cell modality eval. Same recipe as eval_oof_metrics.py
(per-modality StandardScaler -> diff_select(500) -> LinearSVC -> percentile -> donor-grouped 3-fold
CV-OOF -> ridge-NNLS fusion), but adds, per mutation:

  (1) NESTED threshold  — F1-max threshold chosen on the OTHER donor folds, applied to the held fold,
                          so the operating point never sees the samples it scores (kills the optimism).
  (2) BOOTSTRAP 95% CIs — donor-level resampling (B=1000) of AUROC and of sens/spec/F1.
  (3) PERMUTATION null   — per mutation, label shuffles (P=1000) of the fused OOF score's AUROC -> p.
  (4) MODALITY LADDER    — arms: bulk RNA / RNA+Composition / measured-only (drop the 4 imputed-from-RNA
                          blocks ADT,Lipid,Metabolite,GRN) / all-8. Isolates how much multimodal gain is
                          independent of RNA-imputation.

Writes scratchpad/oof_metrics_v2.json.
  bsub -q test -W 150 -M 32000 -R "rusage[mem=32000]" -o ev2.log \
    /usr/local/anaconda3-2020/bin/python eval_oof_metrics_v2.py
"""
import os, sys, json, warnings, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from amlmm.context import build_context, Config
from amlmm import discovery as D, dataio, genetics, udon_features as UF
from amlmm.predictor import diff_select, _pct

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "scratchpad", "oof_metrics_v2.json")
SCORES_OUT = os.path.join(HERE, "..", "scratchpad", "oof_scores_v2.json")   # per-sample (score,y,donor) for calibration + DCA
os.makedirs(os.path.dirname(OUT), exist_ok=True)
SCORES = {}
SC8 = ["RNA", "Composition", "ADT", "Lipid", "Metabolite", "GRN", "LSC", "Cell-comm"]
IMPUTED = {"ADT", "Lipid", "Metabolite", "GRN"}                 # the RNA-conditioned blocks (marked ᴿ in the UI)
MEASURED = [m for m in SC8 if m not in IMPUTED]                 # RNA, Composition, LSC, Cell-comm
ARMS = {"bulkrna": ["BulkRNA"], "rna_comp": ["RNA", "Composition"],
        "measured": MEASURED, "multimodal": SC8}
B_BOOT, P_PERM = 1000, 1000

ctx = build_context(Config(run_id="single_modality"))
samples = ctx.tables["samples"]; dg = samples["donor_group"].astype(str); hold = set(ctx.holdout)
def log(m): print(m, flush=True)


def load_block(mod):
    if mod == "RNA":
        return UF.canonical_rna(ctx)
    c = ctx.path("_sl_%s.pkl" % mod)
    if os.path.exists(c):
        return pd.read_pickle(c)
    if mod == "BulkRNA":
        import bulk_features as BF
        b = BF.bulk_rna_matrix(ctx); b.to_pickle(c); return b
    if mod == "Composition":
        return D._sample_level_matrix(ctx, "composition", set(samples.index))
    if mod in ("ADT", "GRN"):
        return dataio.sample_modality_matrix(ctx, mod)
    if mod in ("Lipid", "Metabolite"):
        return dataio.sample_modality_matrix(ctx, mod, min_spearman=0.3)
    if mod == "Cell-comm":
        return dataio.cellcomm_matrix(ctx)
    if mod == "LSC":
        t = ctx.tables.get("lsc_calls")
        cols = [c for c in ["Prob_m-LSC", "Prob_p+m-LSC", "Prob_p-LSC", "MaxProb"] if c in t.columns]
        return t[cols].apply(pd.to_numeric, errors="coerce")


BLK = {}
for m in sorted(set(sum(ARMS.values(), []))):
    try:
        b = load_block(m).fillna(0.0); BLK[m] = b[~b.index.duplicated(keep="first")]
        log("loaded block %-12s %s" % (m, b.shape))
    except Exception as e:
        log("skip %s: %s" % (m, e))

M = ctx.tables.get("mutations") or genetics.build_mutation_matrix(ctx)
_m01 = {"present": 1.0, "absent": 0.0}


def auc(y, s):
    s = np.asarray(s, float); ok = ~np.isnan(s)
    return roc_auc_score(y[ok], s[ok]) if (ok.sum() >= 4 and len(set(y[ok])) == 2) else np.nan


allidx = set(BLK["RNA"].index) if "RNA" in BLK else set().union(*[set(BLK[m].index) for m in BLK])
MUTS = []
for f in sorted(c for c in M.columns if str(c).startswith(("mut_", "cyto_"))):
    ym = D.labels_for_field(ctx, f).map(_m01)
    tr = [s for s in allidx if pd.notna(ym.get(s)) and s not in hold]
    yv = np.array([int(ym[s]) for s in tr])
    if (yv == 1).sum() >= 8 and (yv == 0).sum() >= 8:
        MUTS.append(f)
log("trainable mutations: %d" % len(MUTS))


def cv_oof(B, ids, y, grp):
    X = B.loc[ids].values; keep = X.std(0) > 0
    if keep.sum() < 2 or len(set(grp)) < 2:
        return None
    oof = np.full(len(ids), np.nan)
    for tri, vai in GroupKFold(min(3, len(set(grp)))).split(X, y, grp):
        if len(set(y[tri])) < 2:
            continue
        sc = StandardScaler().fit(X[tri][:, keep])
        Ztr = sc.transform(X[tri][:, keep]); Zva = sc.transform(X[vai][:, keep])
        sel = diff_select(Ztr, y[tri], 500)
        d = LinearSVC(C=0.02, class_weight="balanced", max_iter=4000).fit(
            Ztr[:, sel], y[tri]).decision_function(Zva[:, sel])
        oof[vai] = d
    ok = ~np.isnan(oof)
    if ok.sum() < 4:
        return None
    p = np.full(len(ids), np.nan); p[ok] = _pct(np.sort(oof[ok]), oof[ok])
    a = auc(y[ok], p[ok])
    if a == a and a < 0.5:
        p = 1 - p
    return p, (max(a, 1 - a) if a == a else 0.5)


def conf(y, call):
    y = np.asarray(y).astype(int); call = np.asarray(call).astype(int)
    tp = int(((call == 1) & (y == 1)).sum()); fp = int(((call == 1) & (y == 0)).sum())
    fn = int(((call == 0) & (y == 1)).sum()); tn = int(((call == 0) & (y == 0)).sum())
    se = tp/(tp+fn) if tp+fn else 0.0; sp = tn/(tn+fp) if tn+fp else 0.0
    pr = tp/(tp+fp) if tp+fp else 0.0; f1 = 2*pr*se/(pr+se) if pr+se else 0.0
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, sensitivity=round(se,4), specificity=round(sp,4),
                precision=round(pr,4), f1=round(f1,4))


def f1max_thr(comb, y):
    bt, bf = 0.5, -1.0
    for t in np.unique(comb):
        c = conf(y, comb >= t)
        if c["f1"] > bf:
            bf, bt = c["f1"], float(t)
    return bt


def nested_metrics(comb, y, groups):
    """threshold chosen on the OTHER folds, applied to the held fold -> optimism-corrected."""
    ng = min(3, len(set(groups)))
    if ng < 2:
        return conf(y, comb >= f1max_thr(comb, y))
    call = np.zeros(len(y), int)
    for tri, tei in GroupKFold(ng).split(comb, y, groups):
        t = f1max_thr(comb[tri], y[tri]) if len(set(y[tri])) == 2 else 0.5
        call[tei] = (comb[tei] >= t).astype(int)
    return conf(y, call)


def boot_ci(comb, y, groups, thr):
    """donor-level bootstrap; AUROC threshold-free, sens/spec/F1 at the fixed operating point."""
    rng = np.random.default_rng(0)
    d2i = {}
    for i, g in enumerate(groups):
        d2i.setdefault(str(g), []).append(i)
    donors = list(d2i); aus, ses, sps, f1s = [], [], [], []
    for _ in range(B_BOOT):
        pick = rng.integers(0, len(donors), len(donors))
        idx = np.concatenate([d2i[donors[k]] for k in pick])
        yb = y[idx]; cb = comb[idx]
        if len(set(yb)) < 2:
            continue
        try:
            aus.append(roc_auc_score(yb, cb))
        except Exception:
            pass
        c = conf(yb, cb >= thr)
        ses.append(c["sensitivity"]); sps.append(c["specificity"]); f1s.append(c["f1"])
    def q(v):
        return [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)] if len(v) > 20 else None
    return dict(auroc_ci=q(aus), sensitivity_ci=q(ses), specificity_ci=q(sps), f1_ci=q(f1s))


def perm_p(comb, y):
    """permutation null: shuffle labels vs the FIXED fused OOF score -> AUROC p (one-sided)."""
    try:
        obs = roc_auc_score(y, comb)
    except Exception:
        return None, None
    rng = np.random.default_rng(1); null = np.empty(P_PERM)
    for i in range(P_PERM):
        null[i] = roc_auc_score(rng.permutation(y), comb)
    p = (1 + int((null >= obs).sum())) / (1 + P_PERM)
    return round(float(obs), 4), round(float(p), 5)


results = {}
for arm, mods in ARMS.items():
    mods = [m for m in mods if m in BLK]
    log("\n===== ARM %s : %s =====" % (arm, mods))
    per_mut = {}
    for mflag in MUTS:
        mshort = mflag.replace("mut_", "").replace("cyto_", "")
        yall = D._labels_for_field_raw(ctx, mflag).map(_m01)
        ym = D.labels_for_field(ctx, mflag).map(_m01)
        oof_p, cvauc, mod_names = {}, {}, []
        for mod in mods:
            B = BLK[mod]
            tr = [s for s in B.index if pd.notna(ym.get(s)) and s not in hold]
            ytr = np.array([int(yall[s]) for s in tr])
            if (ytr == 1).sum() < 5 or (ytr == 0).sum() < 5:
                continue
            mod_names.append(mod)
            cv = cv_oof(B, tr, ytr, dg.loc[tr].values)
            if cv:
                oof_p[mod] = dict(zip(tr, cv[0])); cvauc[mod] = cv[1]
        if not oof_p:
            continue
        common = sorted(set.intersection(*[set(oof_p[m]) for m in oof_p]))
        yco = np.array([int(yall[s]) for s in common]) if common else np.array([])
        if len(common) < 8 or len(set(yco)) != 2:
            continue
        gco = dg.loc[common].values
        cols = [m for m in mod_names if m in oof_p]
        O = np.column_stack([[oof_p[m].get(s, 0.5) for s in common] for m in cols])
        A = np.vstack([O, np.eye(O.shape[1])]); b = np.concatenate([yco.astype(float), np.zeros(O.shape[1])])
        w, _ = nnls(A, b)
        single = {m: auc(yco, np.array([oof_p[m][s] for s in common])) for m in cols}
        bs = max(single, key=lambda k: (single[k] if single[k] == single[k] else -1))
        comb_chk = auc(yco, O @ w / (w.sum() or 1)) if w.sum() > 0 else np.nan
        if not (comb_chk == comb_chk and comb_chk >= (single[bs] if single[bs] == single[bs] else -1) - 1e-9):
            w = np.zeros(len(cols)); w[cols.index(bs)] = 1.0
        wm = {m: round(float(wi / (w.sum() or 1.0)), 3) for m, wi in zip(cols, w)}
        if sum(wm.values()) == 0:
            info = [m for m in mod_names if cvauc.get(m, 0) > 0.52] or mod_names
            wm = {m: round(1.0/len(info), 3) for m in info}
        wc = [m for m in cols if wm.get(m, 0) > 0]
        if not wc:
            continue
        Oc = np.column_stack([[oof_p[m].get(s, 0.5) for s in common] for m in wc])
        wv = np.array([wm[m] for m in wc]); comb = Oc @ wv / (wv.sum() or 1.0)

        naive_t = f1max_thr(comb, yco)
        naive = conf(yco, comb >= naive_t)
        nested = nested_metrics(comb, yco, gco)
        a = auc(yco, comb)
        try:
            ap = float(average_precision_score(yco, comb))
        except Exception:
            ap = float("nan")
        ci = boot_ci(comb, yco, gco, naive_t)
        pobs, pp = perm_p(comb, yco)
        rec = dict(mutation=mshort, n=len(common), n_pos=int((yco == 1).sum()), n_neg=int((yco == 0).sum()),
                   prevalence=round(float(yco.mean()), 4), threshold=round(naive_t, 4),
                   auroc=(round(float(a), 4) if a == a else None), auprc=(round(ap, 4) if ap == ap else None),
                   auprc_baseline=round(float(yco.mean()), 4),
                   sensitivity=naive["sensitivity"], specificity=naive["specificity"], f1=naive["f1"],
                   precision=naive["precision"], tp=naive["tp"], fp=naive["fp"], fn=naive["fn"], tn=naive["tn"],
                   nested_sensitivity=nested["sensitivity"], nested_specificity=nested["specificity"],
                   nested_f1=nested["f1"], nested_precision=nested["precision"],
                   perm_auroc=pobs, perm_p=pp, weights={k: v for k, v in wm.items() if v > 0})
        rec.update(ci)
        per_mut[mshort] = rec
        SCORES.setdefault(arm, {})[mshort] = {          # per-sample OOF for reliability/Brier/ECE + net-benefit
            "score": [round(float(x), 5) for x in comb], "y": [int(v) for v in yco],
            "donor": [str(g) for g in gco]}
        log("  %-14s n+=%-3d auroc=%.3f%s F1n=%.3f F1nest=%.3f permp=%s"
            % (mshort, rec["n_pos"], a, (" CI[%.2f,%.2f]" % tuple(ci["auroc_ci"])) if ci["auroc_ci"] else "",
               naive["f1"], nested["f1"], pp))
    results[arm] = {"modalities": mods, "mutations": per_mut}

# ---- summaries ----
shared = sorted(set.intersection(*[set(results[a]["mutations"]) for a in ARMS]))
ladder = {}
for arm in ARMS:
    rows = [results[arm]["mutations"][m] for m in shared]
    ladder[arm] = {
        "mean_auroc": round(float(np.mean([r["auroc"] for r in rows if r["auroc"] is not None])), 4),
        "mean_nested_f1": round(float(np.mean([r["nested_f1"] for r in rows])), 4),
        "mean_nested_sensitivity": round(float(np.mean([r["nested_sensitivity"] for r in rows])), 4),
        "mean_nested_specificity": round(float(np.mean([r["nested_specificity"] for r in rows])), 4),
        "n_significant_perm_p<0.05": int(sum(1 for r in rows if r["perm_p"] is not None and r["perm_p"] < 0.05)),
        "n_mutations": len(rows)}

payload = {"generated": time.strftime("%Y-%m-%d %H:%M"), "arms": results,
           "shared_mutations": shared, "modality_ladder": ladder,
           "imputed_modalities_dropped_in_measured": sorted(IMPUTED),
           "methods": {"nested_threshold": "F1-max chosen on other donor folds, applied to held fold",
                       "bootstrap": "donor-level resample x%d; AUROC threshold-free, sens/spec/F1 at fixed operating point" % B_BOOT,
                       "permutation": "label shuffle x%d vs fixed fused OOF AUROC; one-sided p" % P_PERM,
                       "caveat": "fusion weights/orientation held fixed in the permutation (a fully re-fused "
                                 "permutation would be more conservative); nested threshold removes the "
                                 "operating-point optimism but the fusion is still label-aware."}}
json.dump(payload, open(OUT, "w"), indent=1)
json.dump(SCORES, open(SCORES_OUT, "w"))
log("wrote per-sample scores -> %s" % SCORES_OUT)
log("\nMODALITY LADDER (mean over %d shared mutations):" % len(shared))
for arm in ARMS:
    L = ladder[arm]
    log("  %-11s AUROC=%.3f  nested F1=%.3f sens=%.3f spec=%.3f  perm-sig %d/%d"
        % (arm, L["mean_auroc"], L["mean_nested_f1"], L["mean_nested_sensitivity"],
           L["mean_nested_specificity"], L["n_significant_perm_p<0.05"], L["n_mutations"]))
log("wrote %s" % OUT)
log("EVAL V2 OK")
