#!/usr/bin/env python3
"""The validation battery for COMPASS-AML Model A.

The failure mode this exists to catch is a model that looks excellent because the observations are not
independent, or because the metric being reported is not the metric the deployed system needs.

Two different tasks are evaluated, and they are NOT interchangeable:

  * **per-drug, across patients** -- "given this drug, which patients respond?" This is what a Spearman
    or AUROC computed within a drug column measures, and it is the number a biomarker paper reports.
  * **per-patient, across drugs** -- "given this patient, which drug should we try?" This is what the
    deployed system is actually asked, and it is harder: it requires the predictions to be comparable
    ACROSS inhibitors, not just ordered within one.

Reported here: continuous (Spearman, RMSE in raw AUC units), tail classification (AUROC, AUPRC with its
own baseline), calibration (ECE, Brier), per-patient top-k retrieval and ranking concordance,
uncertainty-based abstention (coverage vs error), differentiation-state strata, leave-wave-out,
leave-centre-out, a donor-permutation null, and the approved/actionable subset separately.

  python eval_drug_model.py  ->  deliverables/drug_model_validation.json
"""
import os, sys, json, time, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from scipy.stats import spearmanr

from amlmm.drug import data as D, model as M, features as F, targets as TG
from train_drug_model import build_space, fit_predict, sym2ens_map, BUNDLE

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "drug_model_validation.json")
BLOCKS = ["rna", "state", "mut", "clin"]


# --------------------------------------------------------------- helpers ----
def load_all(tail=0.20):
    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(g) for g in d["genes"]]
    ba = [str(s) for s in d["ba_samples"]]
    X = d["ba_X"].astype(np.float64)
    mut = pd.DataFrame(d["ba_L"].astype(float), index=ba, columns=[str(c) for c in d["drivers"]])
    row_of = {s: i for i, s in enumerate(ba)}
    long = D.load(specimens=ba, tail=tail)
    elig = D.eligible_drugs(D.drug_summary(long))
    drugs = sorted(elig.loc[elig["tier"].isin(["primary", "wave_conditional"]), "inhibitor"])
    long = long[long["inhibitor"].isin(drugs)].reset_index(drop=True)
    meta = long.drop_duplicates("specimen").set_index("specimen").reindex(ba).reset_index()
    return dict(genes=genes, ba=ba, X=X, mut=mut, row_of=row_of, long=long, drugs=drugs,
                meta=meta, elig=elig, sym2ens=sym2ens_map(genes))


def oof_predictions(A, folds=5, blocks=BLOCKS, group_col="subject", seed=0, split=None, tail=0.20):
    """Grouped OOF predictions over the whole cohort (or an explicit train/test split).

    `split` = (train_mask, test_mask) runs a single pass instead of k folds -- used for leave-wave-out
    and leave-centre-out, where the partition is given rather than random.
    """
    long, row_of = A["long"], A["row_of"]
    pred = pd.Series(np.nan, index=long.index, dtype=float)
    ycol = pd.Series(np.nan, index=long.index, dtype=float)
    scol = pd.Series(np.nan, index=long.index, dtype=float)
    if split is not None:
        parts = [(np.where(split[0])[0], np.where(split[1])[0])]
    else:
        parts = list(GroupKFold(n_splits=folds).split(long, groups=long[group_col].values))
    for tri, tei in parts:
        f_tr, f_te = long.iloc[tri], long.iloc[tei]
        rows_tr = sorted({row_of[s] for s in f_tr["specimen"]})
        if len(rows_tr) < 40:
            continue
        fs, Z = build_space(A["X"], rows_tr, A["meta"], A["mut"], 100, 4000, A["genes"])
        mod, ltr, lte, p, _ = fit_predict(f_tr, f_te, Z, row_of, A["drugs"], fs, A["sym2ens"],
                                          blocks, A["meta"], A["mut"], tail)
        pred.iloc[tei] = p.values
        ycol.iloc[tei] = lte["y_sens"].values
        scol.iloc[tei] = lte["sens"].values
    out = long.copy()
    out["pred"], out["y_sens"], out["sens"] = pred, ycol, scol
    return out


def per_drug_metrics(df, subset=None):
    rows = []
    for drug, g in df.groupby("inhibitor"):
        if subset is not None and drug not in subset:
            continue
        m = g["pred"].notna() & g["sens"].notna()
        if int(m.sum()) < 20:
            continue
        y, p = g.loc[m, "sens"].values, g.loc[m, "pred"].values
        r = spearmanr(y, p)
        # RMSE back on the raw AUC scale: sens = -(auc - median)/scale, so auc_hat = median - scale*pred
        auc = pd.to_numeric(g.loc[m, "auc"], errors="coerce").values
        med, sc = np.nanmedian(auc), 1.4826 * np.nanmedian(np.abs(auc - np.nanmedian(auc)))
        rmse_auc = float(np.sqrt(np.nanmean((auc - (med - sc * p)) ** 2))) if sc > 0 else np.nan
        e = {"inhibitor": drug, "n": int(m.sum()), "spearman": float(r.statistic),
             "spearman_p": float(r.pvalue), "rmse_z": float(np.sqrt(np.mean((y - p) ** 2))),
             "rmse_auc": rmse_auc, "auc_sd": float(np.nanstd(auc))}
        k = g["y_sens"].notna() & m
        if int(k.sum()) >= 20 and g.loc[k, "y_sens"].nunique() == 2:
            yy = g.loc[k, "y_sens"].astype(int).values; pp = g.loc[k, "pred"].values
            e.update({"n_tail": int(k.sum()), "n_sensitive": int(yy.sum()),
                      "auroc": float(roc_auc_score(yy, pp)),
                      "auprc": float(average_precision_score(yy, pp)),
                      "auprc_baseline": float(yy.mean())})
        rows.append(e)
    return pd.DataFrame(rows)


def per_patient_metrics(df, topk=(1, 3, 5, 10)):
    """The deployment task: within one specimen, rank the inhibitors.

    `hit@k` asks whether any of the model's top-k drugs is in that specimen's true most-sensitive
    decile. `random@k` is the matched chance rate given how many drugs that specimen was tested on --
    quoted alongside every hit rate, because a hit@10 of 0.6 means nothing without it.
    """
    rows = []
    for spec, g in df.groupby("specimen"):
        m = g["pred"].notna() & g["sens"].notna()
        g = g[m]
        if len(g) < 20:
            continue
        obs, prd = g["sens"].values, g["pred"].values
        r = spearmanr(obs, prd)
        n = len(g)
        true_top = set(np.argsort(obs)[::-1][:max(1, int(round(0.10 * n)))])
        e = {"specimen": spec, "n_drugs": n, "spearman": float(r.statistic)}
        order = np.argsort(prd)[::-1]
        from math import comb
        t = len(true_top)
        for k in topk:
            e["hit@%d" % k] = float(len(true_top & set(order[:k])) > 0)
            # exact chance rate for a random k-subset: 1 - C(n-t, k)/C(n, k)
            e["random@%d" % k] = float(1.0 - comb(n - t, k) / comb(n, k)) if n - t >= k else 1.0
            e["prec@%d" % k] = float(len(true_top & set(order[:k])) / k)
        rows.append(e)
    return pd.DataFrame(rows)


def calibration(df, mod_calib, bins=10):
    """ECE / Brier of the per-drug Platt-calibrated probabilities on the tail task."""
    ps, ys = [], []
    for drug, g in df.groupby("inhibitor"):
        k = g["y_sens"].notna() & g["pred"].notna()
        if not int(k.sum()):
            continue
        # go through the model's own calibrator so the percentile step is applied exactly as deployed
        pv = [mod_calib(drug, float(s)) for s in g.loc[k, "pred"].values]
        if any(v is None for v in pv):
            continue
        ps.append(np.asarray(pv, float))
        ys.append(g.loc[k, "y_sens"].astype(int).values)
    if not ps:
        return {}
    p, y = np.concatenate(ps), np.concatenate(ys)
    edges = np.linspace(0, 1, bins + 1)
    ece, curve = 0.0, []
    for i in range(bins):
        m = (p >= edges[i]) & (p < edges[i + 1] if i < bins - 1 else p <= edges[i + 1])
        if not m.sum():
            continue
        ece += m.mean() * abs(p[m].mean() - y[m].mean())
        curve.append({"bin": round(float((edges[i] + edges[i + 1]) / 2), 3),
                      "n": int(m.sum()), "predicted": round(float(p[m].mean()), 4),
                      "observed": round(float(y[m].mean()), 4)})
    return {"n": int(len(y)), "ece": round(float(ece), 4),
            "brier": round(float(brier_score_loss(y, p)), 4),
            "brier_baseline": round(float(brier_score_loss(y, np.full_like(p, y.mean()))), 4),
            "reliability": curve}


def abstention(df, mod_calib):
    """Coverage vs error when the system is allowed to decline. Confidence = |p - 0.5|."""
    ps, ys = [], []
    for drug, g in df.groupby("inhibitor"):
        k = g["y_sens"].notna() & g["pred"].notna()
        if not int(k.sum()):
            continue
        # go through the model's own calibrator so the percentile step is applied exactly as deployed
        pv = [mod_calib(drug, float(s)) for s in g.loc[k, "pred"].values]
        if any(v is None for v in pv):
            continue
        ps.append(np.asarray(pv, float))
        ys.append(g.loc[k, "y_sens"].astype(int).values)
    if not ps:
        return []
    p, y = np.concatenate(ps), np.concatenate(ys)
    conf = np.abs(p - 0.5)
    out = []
    for cov in (1.0, 0.9, 0.75, 0.5, 0.25, 0.1):
        thr = np.quantile(conf, 1.0 - cov)
        m = conf >= thr
        if m.sum() < 20:
            continue
        out.append({"coverage": cov, "n": int(m.sum()),
                    "error_rate": round(float(np.mean((p[m] >= 0.5).astype(int) != y[m])), 4),
                    "auroc": round(float(roc_auc_score(y[m], p[m])), 4) if len(set(y[m])) == 2 else None})
    return out


def strata(df, A):
    """Does performance survive inside differentiation-state strata? If a model only separates
    monocytic from primitive AML it will look good overall and be useless within a stratum."""
    fs = F.FeatureSpace(A["genes"])
    Zall = fs.add_reference("beataml", A["X"])
    S, names = fs._state_block(Zall)
    axis = pd.Series(S[:, names.index("axis_primitive")] - S[:, names.index("axis_mature")],
                     index=A["ba"])
    q = axis.rank(pct=True)
    lab = pd.cut(q, [0, 1 / 3, 2 / 3, 1.0], labels=["monocytic/mature", "intermediate", "primitive"])
    df = df.copy()
    df["stratum"] = df["specimen"].map(lab.to_dict())
    out = {}
    for st, g in df.groupby("stratum", observed=True):
        pd_ = per_drug_metrics(g)
        pp = per_patient_metrics(g)
        out[str(st)] = {"n_specimens": int(g["specimen"].nunique()), "n_drugs": int(len(pd_)),
                        "mean_spearman": round(float(pd_["spearman"].mean()), 4) if len(pd_) else None,
                        "mean_auroc": round(float(pd_["auroc"].mean()), 4) if "auroc" in pd_ else None,
                        "patient_mean_spearman": round(float(pp["spearman"].mean()), 4) if len(pp) else None}
    return out


def permutation_null(A, folds=3, n_perm=20, seed=0, tail=0.20):
    """Break the expression<->response link and nothing else.

    Each specimen is re-pointed at a DIFFERENT specimen's expression row, the same way for every drug.
    That leaves the response matrix, the drug structure, the subject grouping and the class balance
    exactly as they were, so anything the model still finds is an artefact of the procedure rather than
    of the biology. Shuffling AUC within each drug independently would additionally destroy the
    within-patient correlation across drugs, making the null too easy to beat.
    """
    rng = np.random.RandomState(seed)
    real = oof_predictions(A, folds=folds, tail=tail)
    obs = float(per_drug_metrics(real)["spearman"].mean())
    null = []
    specs = np.array(sorted(A["long"]["specimen"].unique()))
    for _ in range(n_perm):
        B = dict(A)
        shuffled = rng.permutation(specs)
        B["row_of"] = {s: A["row_of"][t] for s, t in zip(specs, shuffled)}
        p = oof_predictions(B, folds=folds, tail=tail)
        null.append(float(per_drug_metrics(p)["spearman"].mean()))
    null = np.array(null)
    sd = float(null.std()) or 1e-9
    return {"observed_mean_spearman": round(obs, 4), "n_perm": int(n_perm),
            "null_mean": round(float(null.mean()), 4), "null_sd": round(float(null.std()), 4),
            "null_max": round(float(null.max()), 4),
            "p_value": float((np.sum(null >= obs) + 1) / (n_perm + 1)),
            "p_value_floor": round(1.0 / (n_perm + 1), 4),
            "z_vs_null": round((obs - float(null.mean())) / sd, 2),
            "note": ("with n_perm permutations the smallest attainable p is 1/(n_perm+1); the "
                     "z-distance from the null is the informative statistic when the effect is far "
                     "outside the null range")}


# ------------------------------------------------------------------ main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--perm", type=int, default=20)
    ap.add_argument("--skip-perm", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    A = load_all()
    print("cohort: %d drugs, %d specimens, %d subjects, %d measurements"
          % (len(A["drugs"]), A["long"]["specimen"].nunique(),
             A["long"]["subject"].nunique(), len(A["long"])))

    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "n_drugs": len(A["drugs"]),
           "n_specimens": int(A["long"]["specimen"].nunique()),
           "n_subjects": int(A["long"]["subject"].nunique()),
           "n_measurements": int(len(A["long"]))}

    # ---- main donor-grouped CV over the whole cohort ----
    df = oof_predictions(A, folds=a.folds)
    print("  donor-grouped %d-fold CV done (%.0fs)" % (a.folds, time.time() - t0))
    pd_all = per_drug_metrics(df)
    pp_all = per_patient_metrics(df)

    import pickle
    with open(os.path.join(HERE, "drug_response_model.pkl"), "rb") as f:
        mod = pickle.load(f)

    res["per_drug"] = {
        "n_drugs_scored": int(len(pd_all)),
        "mean_spearman": round(float(pd_all["spearman"].mean()), 4),
        "median_spearman": round(float(pd_all["spearman"].median()), 4),
        "frac_spearman_p_lt_0.05": round(float((pd_all["spearman_p"] < 0.05).mean()), 4),
        "mean_auroc": round(float(pd_all["auroc"].mean()), 4),
        "mean_auprc": round(float(pd_all["auprc"].mean()), 4),
        "mean_auprc_baseline": round(float(pd_all["auprc_baseline"].mean()), 4),
        "mean_rmse_auc": round(float(pd_all["rmse_auc"].mean()), 4),
        "mean_auc_sd": round(float(pd_all["auc_sd"].mean()), 4),
        "table": pd_all.round(4).to_dict("records"),
    }
    res["per_patient"] = {
        "n_specimens": int(len(pp_all)),
        "mean_spearman": round(float(pp_all["spearman"].mean()), 4),
        "median_spearman": round(float(pp_all["spearman"].median()), 4),
        **{k: round(float(pp_all[k].mean()), 4) for k in pp_all.columns
           if k.startswith(("hit@", "random@", "prec@"))},
    }
    res["calibration"] = calibration(df, lambda d, v: mod.calibrated(d, v, "beataml"))
    res["abstention"] = abstention(df, lambda d, v: mod.calibrated(d, v, "beataml"))

    # ---- clinically actionable subset ----
    act = {d for d in A["drugs"] if TG.get(d)["clinical_tier"] in ("approved_AML", "approved_other")}
    pda = per_drug_metrics(df, subset=act)
    res["actionable_subset"] = {
        "n_drugs": int(len(pda)), "definition": "clinical_tier in {approved_AML, approved_other}",
        "mean_spearman": round(float(pda["spearman"].mean()), 4),
        "mean_auroc": round(float(pda["auroc"].mean()), 4),
        "mean_auprc": round(float(pda["auprc"].mean()), 4),
        "table": pda.round(4).sort_values("auroc", ascending=False).to_dict("records"),
    }
    aml = {d for d in A["drugs"] if TG.get(d)["clinical_tier"] == "approved_AML"}
    pdm = per_drug_metrics(df, subset=aml)
    res["approved_AML_subset"] = {"n_drugs": int(len(pdm)),
                                  "table": pdm.round(4).sort_values("auroc", ascending=False).to_dict("records")}

    # ---- differentiation-state strata ----
    res["strata"] = strata(df, A)
    print("  strata done (%.0fs)" % (time.time() - t0))

    # ---- leave-wave-out ----
    coh = A["long"]["cohort"].astype(str)
    w12 = coh.str.contains("1", na=False); w34 = coh.str.contains("3", na=False)
    res["leave_wave_out"] = {}
    for name, tr, te in [("train_w12_test_w34", w12 & ~w34, w34 & ~w12),
                         ("train_w34_test_w12", w34 & ~w12, w12 & ~w34)]:
        if tr.sum() < 2000 or te.sum() < 500:
            continue
        dfw = oof_predictions(A, split=(tr.values, te.values))
        pw = per_drug_metrics(dfw); qw = per_patient_metrics(dfw)
        res["leave_wave_out"][name] = {
            "n_train_rows": int(tr.sum()), "n_test_rows": int(te.sum()),
            "n_drugs": int(len(pw)),
            "mean_spearman": round(float(pw["spearman"].mean()), 4),
            "mean_auroc": round(float(pw["auroc"].mean()), 4) if "auroc" in pw else None,
            "patient_mean_spearman": round(float(qw["spearman"].mean()), 4) if len(qw) else None}
    print("  leave-wave-out done (%.0fs)" % (time.time() - t0))

    # ---- leave-centre-out ----
    ctr = A["long"]["centerID"].astype(str)
    res["leave_center_out"] = {}
    for c, n in ctr.value_counts().head(4).items():
        te = (ctr == c); tr = ~te
        if te.sum() < 500 or tr.sum() < 2000:
            continue
        dfc = oof_predictions(A, split=(tr.values, te.values))
        pc = per_drug_metrics(dfc); qc = per_patient_metrics(dfc)
        res["leave_center_out"]["center_%s" % c] = {
            "n_test_rows": int(te.sum()), "n_drugs": int(len(pc)),
            "mean_spearman": round(float(pc["spearman"].mean()), 4),
            "mean_auroc": round(float(pc["auroc"].mean()), 4) if "auroc" in pc else None,
            "patient_mean_spearman": round(float(qc["spearman"].mean()), 4) if len(qc) else None}
    print("  leave-centre-out done (%.0fs)" % (time.time() - t0))

    # ---- permutation null ----
    if not a.skip_perm:
        res["permutation_null"] = permutation_null(A, folds=3, n_perm=a.perm)
        print("  permutation null done (%.0fs)" % (time.time() - t0))
    elif os.path.exists(OUT):
        prev = json.load(open(OUT)).get("permutation_null")     # 20 min to compute; do not discard it
        if prev:
            res["permutation_null"] = {**prev, "carried_forward_from": "previous run"}

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s  (%.0fs)" % (OUT, time.time() - t0))
    print("\n== headline ==")
    print("  per-drug     mean Spearman %.3f | mean AUROC %.3f | mean AUPRC %.3f (baseline %.3f)"
          % (res["per_drug"]["mean_spearman"], res["per_drug"]["mean_auroc"],
             res["per_drug"]["mean_auprc"], res["per_drug"]["mean_auprc_baseline"]))
    print("  per-patient  mean Spearman %.3f | hit@1 %.3f (chance %.3f) | hit@5 %.3f (chance %.3f)"
          % (res["per_patient"]["mean_spearman"], res["per_patient"]["hit@1"],
             res["per_patient"]["random@1"], res["per_patient"]["hit@5"], res["per_patient"]["random@5"]))
    print("  calibration  ECE %.3f | Brier %.3f (baseline %.3f)"
          % (res["calibration"]["ece"], res["calibration"]["brier"], res["calibration"]["brier_baseline"]))
    print("  actionable   %d drugs, mean AUROC %.3f"
          % (res["actionable_subset"]["n_drugs"], res["actionable_subset"]["mean_auroc"]))


if __name__ == "__main__":
    main()
