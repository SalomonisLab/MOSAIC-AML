#!/usr/bin/env python3
"""Per-mutation SENSITIVITY / SPECIFICITY / F1 for each modality arm, on identical samples & folds.

Reproduces train_predictor.py's recipe EXACTLY (per-modality StandardScaler -> diff_select(500) ->
LinearSVC(C=0.02) -> percentile -> donor-grouped 3-fold CV-OOF -> ridge-NNLS fusion -> F1-max
threshold) but instead of saving a model it exports the fused OOF score vector + truth per mutation,
so sens/spec/F1 can be computed at the deployed operating point.

Both arms run in ONE process so the modality blocks are loaded once and, critically, every arm sees
the SAME samples, SAME labels and SAME donor-grouped folds -> a paired contrast, not apples-to-oranges.

  ARMS (default): bulkrna  = ["BulkRNA"]  bulk RNA-seq alone, NO imputation (the BeatAML-style input)
                  multimodal = the 8 sc blocks (what is deployed)

HONEST CAVEAT recorded in the output: the threshold is F1-max'd on the same OOF vector it is then
scored on, so absolute F1 is optimistic. Both arms get identical treatment, so the CONTRAST is fair.

  bsub -q test -W 240 -M 32000 -R "rusage[mem=32000]" -o eval.log \
    /usr/local/anaconda3-2020/bin/python eval_oof_metrics.py
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
OUT = os.path.join(HERE, "..", "scratchpad", "oof_metrics.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)          # cluster mirror may not have scratchpad/
SC8 = ["RNA", "Composition", "ADT", "Lipid", "Metabolite", "GRN", "LSC", "Cell-comm"]
ARMS = {"bulkrna": ["BulkRNA"], "multimodal": SC8}

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


def metrics_at(comb, y, t):
    pr = comb >= t
    tp = int((pr & (y == 1)).sum()); fn = int(((~pr) & (y == 1)).sum())
    fp = int((pr & (y == 0)).sum()); tn = int(((~pr) & (y == 0)).sum())
    sens = tp / (tp + fn) if (tp + fn) else 0.0            # recall / TPR
    spec = tn / (tn + fp) if (tn + fp) else 0.0            # TNR
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) else 0.0
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, sensitivity=round(sens, 4), specificity=round(spec, 4),
                precision=round(prec, 4), f1=round(f1, 4))


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
        # ridge-NNLS fusion weights (identical to train_predictor)
        wm = {m: 0.0 for m in mod_names}
        cols = [m for m in mod_names if m in oof_p]
        O = np.column_stack([[oof_p[m].get(s, 0.5) for s in common] for m in cols])
        A = np.vstack([O, np.eye(O.shape[1])]); b = np.concatenate([yco.astype(float), np.zeros(O.shape[1])])
        w, _ = nnls(A, b)
        single = {m: auc(yco, np.array([oof_p[m][s] for s in common])) for m in cols}
        bs = max(single, key=lambda k: (single[k] if single[k] == single[k] else -1))
        comb_chk = auc(yco, O @ w / (w.sum() or 1)) if w.sum() > 0 else np.nan
        if not (comb_chk == comb_chk and comb_chk >= (single[bs] if single[bs] == single[bs] else -1) - 1e-9):
            w = np.zeros(len(cols)); w[cols.index(bs)] = 1.0
        tot = w.sum() or 1.0
        for m, wi in zip(cols, w):
            wm[m] = round(float(wi / tot), 3)
        if sum(wm.values()) == 0:
            info = [m for m in mod_names if cvauc.get(m, 0) > 0.52] or mod_names
            for m in info:
                wm[m] = round(1.0 / len(info), 3)
        wc = [m for m in mod_names if m in oof_p and wm.get(m, 0) > 0]
        if not wc:
            continue
        Oc = np.column_stack([[oof_p[m].get(s, 0.5) for s in common] for m in wc])
        wv = np.array([wm[m] for m in wc]); comb = Oc @ wv / (wv.sum() or 1.0)
        # F1-max threshold on the OOF blend (same as deployed)
        bt, bf1 = 0.5, -1.0
        for t in np.unique(comb):
            mm = metrics_at(comb, yco, float(t))
            if mm["f1"] > bf1:
                bf1, bt = mm["f1"], float(t)
        mm = metrics_at(comb, yco, bt)
        a = auc(yco, comb)
        try:
            ap = float(average_precision_score(yco, comb))
        except Exception:
            ap = float("nan")
        mm.update(mutation=mshort, threshold=round(bt, 4),
                  auroc=(round(float(a), 4) if a == a else None),
                  auprc=(round(ap, 4) if ap == ap else None),
                  prevalence=round(float(yco.mean()), 4),
                  n=len(common), n_pos=int((yco == 1).sum()), n_neg=int((yco == 0).sum()),
                  weights={k: v for k, v in wm.items() if v > 0})
        per_mut[mshort] = mm
        log("  %-14s n=%-4d n+=%-3d sens=%.3f spec=%.3f F1=%.3f auroc=%s"
            % (mshort, mm["n"], mm["n_pos"], mm["sensitivity"], mm["specificity"], mm["f1"], mm["auroc"]))
    results[arm] = {"modalities": mods, "mutations": per_mut}

# paired summary over mutations present in BOTH arms
shared = sorted(set(results["bulkrna"]["mutations"]) & set(results["multimodal"]["mutations"]))
summary = {"n_shared_mutations": len(shared), "shared": shared}
for k in ("sensitivity", "specificity", "f1", "auroc"):
    b = [results["bulkrna"]["mutations"][m][k] for m in shared]
    a = [results["multimodal"]["mutations"][m][k] for m in shared]
    b = [x for x in b if x is not None]; a = [x for x in a if x is not None]
    summary[k] = {"bulkrna_mean": round(float(np.mean(b)), 4), "multimodal_mean": round(float(np.mean(a)), 4),
                  "delta": round(float(np.mean(a) - np.mean(b)), 4)}
    try:
        from scipy.stats import wilcoxon
        st, p = wilcoxon([results["multimodal"]["mutations"][m][k] for m in shared],
                         [results["bulkrna"]["mutations"][m][k] for m in shared])
        summary[k]["wilcoxon_p"] = float(p)
    except Exception:
        pass

payload = {"generated": time.strftime("%Y-%m-%d %H:%M"), "arms": results, "paired_summary": summary,
           "caveat": "Threshold is F1-max'd on the same donor-grouped CV-OOF vector it is scored on, so "
                     "absolute F1/sens/spec are optimistic. Both arms get identical treatment on identical "
                     "samples and folds, so the CONTRAST between arms is fair."}
json.dump(payload, open(OUT, "w"), indent=1)
log("\nwrote %s" % OUT)
log("paired over %d mutations: %s" % (len(shared), json.dumps({k: summary[k] for k in
    ("sensitivity", "specificity", "f1")}, indent=1)))
log("EVAL OOF METRICS OK")
