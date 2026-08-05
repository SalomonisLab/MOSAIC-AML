#!/usr/bin/env python3
"""Bulk / bulk+sc mutation-prediction bake-off across cohorts and model families.

4 tests  = {BeatAML, BeatAML+Leucegene} x {bulk-only, bulk+sc}   (bulk+sc pools our scRNA-as-bulk into the
           training set; different patients, shared gene space)
x FS      = {feature-selection ON (diff-select top-K), OFF (all genes)}
= 8 runs, each running the full model panel:
   logL2 . linSVM . shrLDA . RF . HistGB . PLS . MLP(neural net) . ENSEMBLE

EVALUATION (per the request): every run scored on our sealed 29 (26 w/ bulk) single-cell held-out; the
BeatAML-only runs (which never train on Leucegene) ALSO externally validated on Leucegene.

Output: gui/bulk_bakeoff_results.json  (tidy: test, fs, driver, model, eval_target -> AUROC + n).
Env: AMLMM_SMOKE=1 -> few drivers / FS-on only / fast models (validation).  Runs locally.
"""
import os, sys, json, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
from sklearn.utils.class_weight import compute_sample_weight
from scipy.stats import rankdata
from bulk_external import harmonize, DRIVERS
from amlmm.predictor import diff_select

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SMOKE = os.environ.get("AMLMM_SMOKE", "") == "1"
VAR_CAP = int(os.environ.get("AMLMM_VARCAP", "2500"))   # FS-OFF caps to top-N most-variable genes (label-agnostic)
MARK_TOPK = int(os.environ.get("AMLMM_MARK_TOPK", "25"))  # FS-ON markers kept per dataset (10-50 range)
MARK_MAX = int(os.environ.get("AMLMM_MARK_MAX", "50"))    # cap on the per-dataset-union marker panel
def log(m): print(m, flush=True)


def mwu_p_greater(mut, non):
    """vectorized one-sided Mann-Whitney (H1: mut > non) p per gene; normal approx, no scipy-axis dependency."""
    from scipy.stats import norm
    n1, n0 = mut.shape[0], non.shape[0]
    ranks = np.vstack([mut, non]).argsort(0).argsort(0).astype(float) + 1.0
    U1 = ranks[:n1].sum(0) - n1 * (n1 + 1) / 2.0
    sd = np.sqrt(n1 * n0 * (n1 + n0 + 1) / 12.0) or 1.0
    return 1.0 - norm.cdf((U1 - n1 * n0 / 2.0) / sd)


def marker_select(Xlin, y, top_k=MARK_TOPK, fold_thr=2.0, p_thr=0.05):
    """MarkerFinder-style: genes with linear fold>2 AND MWU p<0.05 (mutant vs non-mutant), top-K by fold.
    Xlin = LINEAR expression for one dataset; returns column indices (into the shared common-gene order)."""
    mut, non = Xlin[y == 1], Xlin[y == 0]
    if len(mut) < 3 or len(non) < 3:
        return set()
    eps = 1e-6
    fold = (mut.mean(0) + eps) / (non.mean(0) + eps)
    p = mwu_p_greater(mut, non)
    passing = np.where((fold > fold_thr) & (p < p_thr))[0]
    if not len(passing):
        return set()
    return set(int(i) for i in passing[np.argsort(fold[passing])[::-1]][:top_k])


def score_models(Xtr, ytr, evalXs, fast=False):
    """fit each model on (Xtr,ytr); return {model: {eval_name: score_vector}}."""
    sw = compute_sample_weight("balanced", ytr)
    defs = [
        ("logL2", lambda: LogisticRegression(C=0.05, class_weight="balanced", max_iter=2000), "proba"),
        ("linSVM", lambda: LinearSVC(C=0.02, class_weight="balanced", max_iter=4000), "dec"),
        ("shrLDA", lambda: LinearDiscriminantAnalysis(solver="lsqr", shrinkage=0.5), "proba"),
        ("PLS", lambda: PLSRegression(n_components=min(10, Xtr.shape[1])), "pls"),
        ("RF", lambda: RandomForestClassifier(n_estimators=300, class_weight="balanced_subsample", n_jobs=int(os.environ.get("AMLMM_NJOBS", "4")), random_state=0), "proba"),
        ("HistGB", lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=150, l2_regularization=1.0, random_state=0), "proba_sw"),
        ("MLP", lambda: MLPClassifier(hidden_layer_sizes=(128, 32), alpha=1e-3, early_stopping=True, max_iter=300, random_state=0), "proba"),
    ]
    if fast:
        defs = [d for d in defs if d[0] in ("logL2", "linSVM", "shrLDA", "PLS")]
    out = {}
    for name, mk, kind in defs:
        try:
            est = mk()
            if kind == "proba_sw":
                est.fit(Xtr, ytr, sample_weight=sw); f = lambda Xe: est.predict_proba(Xe)[:, 1]
            elif kind == "pls":
                est.fit(Xtr, ytr.astype(float)); f = lambda Xe: est.predict(Xe).ravel()
            elif kind == "dec":
                est.fit(Xtr, ytr); f = lambda Xe: est.decision_function(Xe)
            else:
                est.fit(Xtr, ytr); f = lambda Xe: est.predict_proba(Xe)[:, 1]
            out[name] = {en: f(Xe) for en, Xe in evalXs.items()}
        except Exception as e:
            out[name] = {en: None for en in evalXs}
    # ensemble = rank-average of the linear members
    ens_mem = [m for m in ("logL2", "linSVM", "shrLDA", "PLS") if m in out]
    out["ENSEMBLE"] = {}
    for en in evalXs:
        rs = [rankdata(out[m][en]) for m in ens_mem if out[m][en] is not None]
        out["ENSEMBLE"][en] = np.mean(rs, axis=0) if rs else None
    return out


def main():
    t0 = time.time()
    log("harmonizing cohorts ...")
    H = harmonize()
    genes = H["genes"]; hold = H["holdout"]
    baZ, baLin, baL = H["beataml"]["Xz"], H["beataml"]["Xlin"], H["beataml"]["L"]
    lgZ, lgLin, lgL = H["leucegene"]["Xz"], H["leucegene"]["Xlin"], H["leucegene"]["L"]
    scZ, scLin, scL = H["sc"]["Xz"], H["sc"]["Xlin"], H["sc"]["L"]
    sc_hold = [s for s in scZ.index if s in hold]
    sc_train = [s for s in scZ.index if s not in hold]
    log("genes %d | BeatAML %d | Leucegene %d | sc-train %d | sc-holdout %d (%.1fs)"
        % (len(genes), len(baZ), len(lgZ), len(sc_train), len(sc_hold), time.time() - t0))

    drivers = ["NPM1", "TET2", "TP53", "FLT3-ITD"] if SMOKE else DRIVERS
    fs_opts = [("fs_on", "markers")] if SMOKE else [("fs_on", "markers"), ("fs_off", None)]
    tests = {"beataml": [(baZ, baLin, baL)],
             "beataml_leucegene": [(baZ, baLin, baL), (lgZ, lgLin, lgL)]}
    rows = []
    outp = os.path.join(ROOT, "gui", "bulk_bakeoff_results.json")
    def dump(done):                                                    # incremental: survives a walltime kill
        json.dump({"rows": rows, "drivers": drivers, "complete": done, "built": time.strftime("%Y-%m-%d %H:%M"),
                   "note": "AUROC by test x modality x fs x driver x model x eval-target"}, open(outp, "w"), indent=1)
    for test, ext_parts in tests.items():
        for mod in ["bulk", "bulk_sc"]:
            parts = list(ext_parts) + ([(scZ.loc[sc_train], scLin.loc[sc_train], scL.loc[sc_train])] if mod == "bulk_sc" else [])
            for fs_name, K in fs_opts:
                for drv in drivers:
                    Zs, ys, per_lin = [], [], []                            # Z = z-scored (classifier); Lin = linear (FS fold)
                    for Z, Lin, L in parts:
                        y = L[drv]; keep = y.notna(); idx = keep.index[keep]
                        if len(idx):
                            yv = y.loc[idx].astype(int).values
                            Zs.append(Z.loc[idx]); ys.append(yv); per_lin.append((Lin.loc[idx].values, yv))
                    if not Zs:
                        continue
                    Ztr_df = pd.concat(Zs); ytr = np.concatenate(ys)
                    if (ytr == 1).sum() < 8 or (ytr == 0).sum() < 8:
                        continue
                    evalXs, evalYs = {}, {}
                    yh = scL.loc[sc_hold, drv]
                    if yh.notna().sum() >= 4 and yh.dropna().nunique() == 2:
                        ix = yh.index[yh.notna()]; evalXs["sc_heldout"] = scZ.loc[ix]; evalYs["sc_heldout"] = yh.loc[ix].astype(int).values
                    if test == "beataml" and drv in lgL:
                        yl = lgL[drv]
                        if yl.notna().sum() >= 8 and yl.dropna().nunique() == 2:
                            ix = yl.index[yl.notna()]; evalXs["leucegene"] = lgZ.loc[ix]; evalYs["leucegene"] = yl.loc[ix].astype(int).values
                    if not evalXs:
                        continue
                    if K == "markers":                                     # FS ON: fold>2 & MWU p<0.05 markers, per-dataset union
                        gset = set()
                        for Xlin_d, y_d in per_lin:
                            gset |= marker_select(Xlin_d, y_d)
                        gsel = sorted(gset)[:MARK_MAX]
                        if len(gsel) < 3:
                            continue
                    else:                                                  # FS OFF: top-VAR_CAP most-variable (tractability only)
                        gsel = list(np.argsort(Ztr_df.values.var(0))[::-1][:VAR_CAP])
                    Ztr = Ztr_df.values[:, gsel]
                    ev = {en: Xe.values[:, gsel] for en, Xe in evalXs.items()}
                    res = score_models(Ztr, ytr, ev, fast=SMOKE)
                    nmark = len(gsel)
                    for model, sc_out in res.items():
                        for en, sv in sc_out.items():
                            au = roc_auc_score(evalYs[en], sv) if sv is not None and len(set(evalYs[en])) == 2 else None
                            rows.append({"test": test, "modality": mod, "fs": fs_name, "driver": drv, "model": model,
                                         "eval": en, "auroc": (round(float(au), 3) if au is not None else None),
                                         "n_feat": int(nmark), "n_train": int(len(ytr)), "n_train_pos": int((ytr == 1).sum()),
                                         "n_eval_pos": int((evalYs[en] == 1).sum())})
                log("done %-18s %-8s %-6s (%.0fs, %d rows)" % (test, mod, fs_name, time.time() - t0, len(rows)))
                dump(False)

    dump(True)
    log("wrote %s (%d rows, %.0fs)" % (outp, len(rows), time.time() - t0))
    # quick summary: mean AUROC on sc_heldout by test/modality/fs/model
    df = pd.DataFrame([r for r in rows if r["eval"] == "sc_heldout" and r["auroc"] is not None])
    if len(df):
        piv = df.groupby(["test", "modality", "fs", "model"])["auroc"].mean().round(3)
        log("\n=== mean held-out AUROC (sc) by config x model ===\n" + piv.to_string())
    log("BULK BAKEOFF OK")


if __name__ == "__main__":
    main()
