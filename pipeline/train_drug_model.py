#!/usr/bin/env python3
"""Train and honestly evaluate COMPASS-AML Model A on the BeatAML2 ex-vivo inhibitor screen.

Split discipline, in order of how easy each is to get wrong:

  1. **Patients, never specimens.** 34 BeatAML subjects contributed >1 drug-tested specimen. Splitting
     on specimens would put a patient's diagnosis sample in train and their relapse sample in test.
  2. **A sealed hold-out of >=15% of PATIENTS**, drawn once with a fixed seed and never touched by any
     fitting, feature selection, threshold or calibration step.
  3. **Everything refits inside every fold** -- the cohort z-reference, the PCA, the clinical
     imputation medians, the within-drug normalisation constants and the sensitive/resistant tail
     cut-points. Normalising AUC against the whole cohort before splitting is the subtle leak that
     makes ex-vivo response models look better than they are.

  python train_drug_model.py [--folds 5] [--holdout 0.15] [--blocks rna,state,mut,clin]
      -> pipeline/drug_response_model.pkl
      -> deliverables/drug_model_card.json  (per-drug OOF + sealed-hold-out metrics)
"""
import os, sys, json, time, pickle, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

from amlmm.drug import data as D, model as M, features as F, targets as TG

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
NORM_EXP = os.path.join(ROOT, "data", "external", "beataml", "norm_exp.txt")
OUTP = os.path.join(HERE, "drug_response_model.pkl")
CARD = os.path.join(ROOT, "deliverables", "drug_model_card.json")


def sym2ens_map(genes):
    hdr = pd.read_csv(NORM_EXP, sep="\t", usecols=["stable_id", "display_label"])
    s = {}
    gs = set(genes)
    for e, y in zip(hdr["stable_id"].astype(str), hdr["display_label"].astype(str)):
        if e in gs:
            s.setdefault(y, e)
    return s


def build_space(X_lin_all, train_rows, meta_all, mut_all, n_pc, var_genes, genes):
    """Fit a FeatureSpace on `train_rows` only, then transform every specimen with it."""
    fs = F.FeatureSpace(genes, n_pc=n_pc, var_genes=var_genes)
    fs.add_reference("beataml", X_lin_all[train_rows])          # reference = TRAIN specimens only
    Z = fs.z(X_lin_all, "beataml")
    fs.fit(Z[train_rows], meta_all.iloc[train_rows], mut_all.iloc[train_rows])
    return fs, Z


def fit_predict(long_tr, long_te, Z, row_of, drugs, fs, sym2ens, blocks, meta_all, mut_all, tail):
    """One (train -> test) pass: normalise on train, build features, fit, predict test rows."""
    norm = M.fit_norm(long_tr, np.ones(len(long_tr), bool), tail=tail)
    ltr = M.apply_norm(long_tr, norm)
    lte = M.apply_norm(long_te, norm)
    Xp, _, sl = fs.transform(Z, meta=meta_all, mut=mut_all, blocks=blocks)
    mod = M.DrugResponseModel(fs).prepare_targets(Z, sym2ens)
    mod.norm = norm
    mod.sym2ens = sym2ens
    mod.fit(ltr, Xp, row_of, drugs, sl)
    return mod, ltr, lte, mod.predict_rows(lte, Xp, row_of), Xp


def metrics(sub, pred):
    """Continuous + tail-classification metrics for one drug."""
    m = pred.notna()
    y = sub.loc[m, "sens"].astype(float)
    p = pred[m].astype(float)
    out = {"n": int(m.sum())}
    if len(y) >= 10 and y.std() > 0:
        r = spearmanr(y, p)
        out["spearman"] = round(float(r.statistic), 4)
        out["spearman_p"] = float(r.pvalue)
        out["rmse_z"] = round(float(np.sqrt(np.mean((y - p) ** 2))), 4)
    ym = sub.loc[m, "y_sens"]
    k = ym.notna()
    if int(k.sum()) >= 10 and ym[k].nunique() == 2:
        yy = ym[k].astype(int).values
        pp = p[k.values].values
        out["n_tail"] = int(k.sum()); out["n_sensitive"] = int(yy.sum())
        out["auroc"] = round(float(roc_auc_score(yy, pp)), 4)
        out["auprc"] = round(float(average_precision_score(yy, pp)), 4)
        out["auprc_baseline"] = round(float(yy.mean()), 4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--holdout", type=float, default=0.15)
    ap.add_argument("--tail", type=float, default=0.20)
    ap.add_argument("--n-pc", type=int, default=100)
    ap.add_argument("--var-genes", type=int, default=4000)
    ap.add_argument("--blocks", default="rna,state,mut,clin")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    blocks = [b for b in a.blocks.split(",") if b]
    t0 = time.time()

    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(g) for g in d["genes"]]
    ba_samples = [str(s) for s in d["ba_samples"]]
    X_lin = d["ba_X"].astype(np.float64)
    mut_all = pd.DataFrame(d["ba_L"].astype(float), index=ba_samples,
                           columns=[str(c) for c in d["drivers"]])
    row_of = {s: i for i, s in enumerate(ba_samples)}
    sym2ens = sym2ens_map(genes)

    long = D.load(specimens=ba_samples, tail=a.tail)
    summ = D.drug_summary(long)
    elig = D.eligible_drugs(summ)
    drugs = sorted(elig.loc[elig["tier"].isin(["primary", "wave_conditional"]), "inhibitor"])
    long = long[long["inhibitor"].isin(drugs)].reset_index(drop=True)

    meta_all = (long.drop_duplicates("specimen").set_index("specimen")
                .reindex(ba_samples).reset_index())

    subs = np.array(sorted(long["subject"].unique()))
    rng = np.random.RandomState(a.seed)
    perm = rng.permutation(len(subs))
    n_hold = int(round(a.holdout * len(subs)))
    hold_subs = set(subs[perm[:n_hold]])
    train_subs = set(subs[perm[n_hold:]])
    is_hold = long["subject"].isin(hold_subs)

    print("COMPASS-AML Model A  |  blocks=%s" % ",".join(blocks))
    print("  drugs %d | specimens %d | subjects %d | measurements %d"
          % (len(drugs), long["specimen"].nunique(), len(subs), len(long)))
    print("  sealed hold-out: %d subjects (%.1f%%), %d specimens, %d measurements"
          % (len(hold_subs), 100.0 * len(hold_subs) / len(subs),
             long.loc[is_hold, "specimen"].nunique(), int(is_hold.sum())))

    # ---------------- donor-grouped CV inside the training partition ----------------
    tr_long = long[~is_hold].reset_index(drop=True)
    oof = pd.Series(np.nan, index=tr_long.index, dtype=float)
    oof_y = pd.Series(np.nan, index=tr_long.index, dtype=float)
    oof_sens = pd.Series(np.nan, index=tr_long.index, dtype=float)
    gkf = GroupKFold(n_splits=a.folds)
    for k, (tri, tei) in enumerate(gkf.split(tr_long, groups=tr_long["subject"].values)):
        f_tr, f_te = tr_long.iloc[tri], tr_long.iloc[tei]
        rows_tr = sorted({row_of[s] for s in f_tr["specimen"]})
        fs, Z = build_space(X_lin, rows_tr, meta_all, mut_all, a.n_pc, a.var_genes, genes)
        mod, ltr, lte, pred, _ = fit_predict(f_tr, f_te, Z, row_of, drugs, fs, sym2ens,
                                             blocks, meta_all, mut_all, a.tail)
        oof.iloc[tei] = pred.values
        oof_y.iloc[tei] = lte["y_sens"].values
        oof_sens.iloc[tei] = lte["sens"].values
        if not a.quiet:
            print("   fold %d/%d  train %d specimens  test %d rows  (%.0fs)"
                  % (k + 1, a.folds, len(rows_tr), len(f_te), time.time() - t0))
    tr_long["y_sens"], tr_long["sens"] = oof_y, oof_sens

    # ---------------- deployed model: all training subjects ----------------
    rows_tr = sorted({row_of[s] for s in tr_long["specimen"]})
    fs, Z = build_space(X_lin, rows_tr, meta_all, mut_all, a.n_pc, a.var_genes, genes)
    ho_long = long[is_hold].reset_index(drop=True)
    mod, ltr_full, lho, pred_ho, Xp = fit_predict(tr_long, ho_long, Z, row_of, drugs, fs, sym2ens,
                                                  blocks, meta_all, mut_all, a.tail)
    mod.fit_calibration(tr_long, oof)                      # calibrate on OOF, never on the hold-out
    mod.set_score_reference(ltr_full, mod.predict_rows(ltr_full, Xp, row_of))
    mod.attach_neighbours(Z[rows_tr], [ba_samples[i] for i in rows_tr], ltr_full)
    mod.blocks = blocks
    mod.sym2ens = sym2ens
    mod.holdout_subjects = sorted(hold_subs)
    mod.drugs = drugs

    # ---------------- metrics ----------------
    per = {}
    for drug in drugs:
        s_oof = tr_long[tr_long["inhibitor"] == drug]
        s_ho = lho[lho["inhibitor"] == drug]
        per[drug] = {"oof": metrics(s_oof, oof.loc[s_oof.index]),
                     "holdout": metrics(s_ho, pred_ho.loc[s_ho.index]),
                     "annotation": {k: v for k, v in TG.get(drug).items()
                                    if k in ("family", "family_group", "clinical_tier", "exposure",
                                             "targets", "mechanism", "analogue")},
                     "tier": elig.loc[elig["inhibitor"] == drug, "tier"].iloc[0]}

    def agg(key, metric):
        v = [per[d][key].get(metric) for d in drugs if per[d][key].get(metric) is not None]
        return (round(float(np.mean(v)), 4), len(v)) if v else (None, 0)

    card = {
        "version": M.DrugResponseModel.VERSION,
        "trained": time.strftime("%Y-%m-%d %H:%M"),
        "blocks": blocks, "folds": a.folds, "tail": a.tail, "seed": a.seed,
        "n_drugs": len(drugs), "n_specimens": int(long["specimen"].nunique()),
        "n_subjects": int(len(subs)), "n_measurements": int(len(long)),
        "holdout": {"subjects": len(hold_subs), "specimens": int(long.loc[is_hold, "specimen"].nunique()),
                    "measurements": int(is_hold.sum()), "fraction_subjects": round(len(hold_subs) / len(subs), 3)},
        "qc": {"rows_before": long.attrs.get("n_rows_before_qc"),
               "dropped_by_flag": long.attrs.get("qc_dropped_by_flag")},
        "summary": {
            "oof_mean_spearman": agg("oof", "spearman"), "oof_mean_auroc": agg("oof", "auroc"),
            "oof_mean_auprc": agg("oof", "auprc"),
            "holdout_mean_spearman": agg("holdout", "spearman"),
            "holdout_mean_auroc": agg("holdout", "auroc"),
            "holdout_mean_auprc": agg("holdout", "auprc"),
        },
        "per_drug": per,
        "model": mod.meta,
    }
    # the deployed model carries its own honest performance, so downstream scoring can down-weight a
    # drug it cannot predict instead of trusting every point estimate equally
    # assay reliability travels WITH the model: a recommendation for an inhibitor whose own
    # measurement does not reproduce must carry that fact to the point of use, not sit in a side file
    rel = D.drug_reliability(long)
    mod.reliability = dict(zip(rel["inhibitor"], rel["reliability"]))
    mod.reliability_tier = dict(zip(rel["inhibitor"], rel["reliability_tier"]))
    card["assay_reliability"] = rel.round(4).to_dict("records")
    card["assay_reliability_summary"] = rel["reliability_tier"].value_counts().to_dict()
    mod.oof_metrics = {d: per[d]["oof"] for d in drugs}
    mod.holdout_metrics = {d: per[d]["holdout"] for d in drugs}
    mod.drug_tier = {d: per[d]["tier"] for d in drugs}

    tagged = CARD if not a.tag else CARD.replace(".json", "_%s.json" % a.tag)
    json.dump(card, open(tagged, "w"), indent=1)
    if not a.tag:
        with open(OUTP, "wb") as f:
            pickle.dump(mod, f)

    print("\n== summary (%d drugs) ==" % len(drugs))
    for k, v in card["summary"].items():
        print("  %-26s %s   (n drugs = %s)" % (k, v[0], v[1]))
    top = sorted(drugs, key=lambda d: -(per[d]["oof"].get("auroc") or 0))[:12]
    print("\n%-26s %6s %6s %6s | %6s %6s   %s" % ("drug", "sp_oof", "auc_oof", "n", "sp_ho", "auc_ho", "tier"))
    for dd in top:
        o, h = per[dd]["oof"], per[dd]["holdout"]
        print("%-26s %6s %6s %6s | %6s %6s   %s"
              % (dd[:26], o.get("spearman"), o.get("auroc"), o.get("n"),
                 h.get("spearman"), h.get("auroc"), per[dd]["annotation"]["clinical_tier"]))
    print("\nwrote %s%s  (%.0fs)" % (tagged, "" if a.tag else " + drug_response_model.pkl", time.time() - t0))


if __name__ == "__main__":
    main()
