#!/usr/bin/env python3
"""Fit a cohort adapter: per mutation, use the cohort-trained model only where it beats the shipped caller.

  python train_cohort_adapter.py --expr <bulk_equiv.npz> --truth <per_patient.tsv> --cohort 32266596 \
      --out cohort_adapters/32266596.json.gz [--exclude "WHO=Healthy control"] [--sealed-frac 0.3]

--expr   npz with X (patients x genes, linear: counts / CP10k / TPM - renormalised internally), genes
         (symbols or ENSG), patients (ids matching --truth). For single-cell data this is the per-patient
         bulk-equivalent, i.e. what ingest_patient.bulk_equiv_from_adata produces.
--truth  per-patient table; mutation columns TRUE/FALSE, blank/NA = not assayed (left UNLABELLED, never
         counted as wild-type).

Honesty guarantees, in order:
  1. A sealed fraction of patients is set aside FIRST (ids recorded in the adapter). Nothing - not the
     gene selection, not the per-unit routing decision, not a threshold - sees them.
  2. The per-unit routing decision (cohort model vs shipped) is made from training folds only.
  3. The whole procedure is scored by nested CV on the development set: each outer fold refits the
     adapter, routing decisions included, and scores patients it never saw. Selecting the winner per
     mutation from the same numbers you report is the optimism this avoids.
  4. The final adapter (fitted on the development set) is scored once on the sealed patients.
"""
from __future__ import annotations
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from amlmm.bulk_predictor import BulkMutationPredictor
from amlmm import cohort_adapter as CA

HERE = os.path.dirname(os.path.abspath(__file__))

# Units for a per-patient table in the 32266596 schema (TRUE/FALSE gene columns, cytogenetic flags).
#   category: replaces one shipped category;  gene: one gene over several variant categories;
#   new: a category the shipped caller does not have.
UNITS_32266596 = (
    [{"name": g, "mode": "category", "replaces": [c], "truth": [t]} for g, c, t in [
        ("NPM1", "NPM1_exon12_frameshift", "NPM1"), ("RUNX1", "RUNX1_LOF", "RUNX1"), ("TET2", "TET2", "TET2"),
        ("ASXL1", "ASXL1", "ASXL1"), ("SRSF2", "SRSF2", "SRSF2"), ("BCOR", "BCOR", "BCOR"),
        ("WT1", "WT1_LOF", "WT1"), ("KIT", "KIT_D816", "KIT"), ("IDH1", "IDH1_R132", "IDH1"),
        ("STAG2", "STAG2", "STAG2"), ("PTPN11", "PTPN11_other", "PTPN11"),
        ("FLT3-ITD", "FLT3_ITD", "FLT3.ITD"), ("FLT3-TKD", "FLT3_TKD_D835/I836", "FLT3.TKD"),
        ("CEBPA-bZIP", "CEBPA_bZIP", "CEBPA_bZIP")]]
    + [{"name": g, "mode": "gene", "replaces": cs, "truth": [g]} for g, cs in [
        ("TP53", ["TP53_hotspot_DBD", "TP53_LOF/splice/frameshift"]),
        ("NRAS", ["NRAS_G12", "NRAS_G13", "NRAS_Q61"]),
        ("DNMT3A", ["DNMT3A_R882", "DNMT3A_nonR882"]),
        ("IDH2", ["IDH2_R140", "IDH2_R172"])]]
    + [{"name": n, "mode": "new", "category": c, "label": lab, "truth": t} for n, c, lab, t in [
        ("t(15;17) APL", "t15_17", "t(15;17) PML::RARA", ["t.15.17"]),
        ("inv(16)", "inv16", "inv(16) CBFB::MYH11", ["inv16"]),
        ("t(8;21)", "t8_21", "t(8;21) RUNX1::RUNX1T1", ["t.8.21"]),
        ("-7", "del7", "monosomy 7", ["Monosomy7"]),
        ("+8", "trisomy8", "trisomy 8", ["Trisomy8"]),
        ("inv(3)", "inv3", "inv(3)/t(3;3) MECOM", ["inv3"]),
        ("complex karyotype", "complex", "complex karyotype", {"col": "broad_group", "pos": ["complex"],
                                                                 "neg": ["normal", "specific"]})]])


def truth_vector(tab, spec):
    if isinstance(spec, dict):
        v = tab[spec["col"]].astype(str)
        return np.where(v.isin(spec["pos"]), 1.0, np.where(v.isin(spec["neg"]), 0.0, np.nan))
    cols = [tab[c].map({"TRUE": 1.0, "FALSE": 0.0, True: 1.0, False: 0.0}) for c in spec]
    M = np.column_stack([c.values.astype(float) for c in cols])
    anypos = np.nanmax(np.where(np.isnan(M), -1, M), axis=1) == 1
    allneg = np.all(M == 0, axis=1)
    return np.where(anypos, 1.0, np.where(allneg, 0.0, np.nan))


def unit_output(u, calls):
    """(score, called) of one unit in one patient's call dict; None if the system has no such output."""
    if u["mode"] == "new":
        r = calls.get(u["category"])
        return (None, None) if r is None else (r["probability"], r["call"] == "present")
    rs = [calls[c] for c in u["replaces"] if c in calls]
    if not rs:
        return None, None
    return max(r["probability"] or 0 for r in rs), any(r["call"] == "present" for r in rs)


def score(units, labels, systems, rows, fold=None):
    """systems: {name: list-of-call-dicts aligned with rows}. Per-unit AUROC + pooled call metrics.

    fold: outer-fold id per row. Scores are rank-normalised WITHIN each fold before pooling: each fold
    has its own adapter, so a unit routed in one fold and not another would otherwise pool two score
    scales (shipped percentile vs adapter percentile) and distort the AUROC. Calls need no such care."""
    from scipy.stats import rankdata
    fold = np.zeros(len(rows), int) if fold is None else np.asarray(fold)
    per, pooled = [], {}
    for s in systems:
        pooled[s] = {"shipped_units": np.zeros(4, int), "new_units": np.zeros(4, int)}
    for u in units:
        y = labels[u["name"]][rows]
        m = np.isfinite(y)
        if m.sum() == 0 or np.nansum(y) == 0:
            continue
        rec = {"unit": u["name"], "mode": u["mode"], "n_pos": int(np.nansum(y)), "n": int(m.sum())}
        for s, calls in systems.items():
            outs = [unit_output(u, calls[i]) for i in range(len(rows))]
            if all(o[0] is None for o in outs):
                rec[s + "_auroc"] = None
                continue
            raw = np.array([o[0] if o[0] is not None else 0.0 for o in outs], float)
            for f in np.unique(fold):
                ff = fold == f
                raw[ff] = rankdata(raw[ff]) / ff.sum()
            sc = raw[m]
            cl = np.array([bool(o[1]) for o in outs])[m]
            yy = y[m].astype(int)
            rec[s + "_auroc"] = round(float(roc_auc_score(yy, sc)), 3) if 0 < yy.sum() < len(yy) else None
            tp, fp = int((cl & (yy == 1)).sum()), int((cl & (yy == 0)).sum())
            fn, tn = int((~cl & (yy == 1)).sum()), int((~cl & (yy == 0)).sum())
            rec[s + "_calls"] = [tp, fp, fn, tn]
            pooled[s]["new_units" if u["mode"] == "new" else "shipped_units"] += [tp, fp, fn, tn]
        per.append(rec)

    def summ(c):
        tp, fp, fn, tn = (int(v) for v in c)
        pr = tp / max(1, tp + fp); se = tp / max(1, tp + fn)
        return {"TP": tp, "FP": fp, "FN": fn, "TN": tn, "sensitivity": round(se, 3),
                "specificity": round(tn / max(1, tn + fp), 3), "precision": round(pr, 3),
                "F1": round(2 * pr * se / max(1e-9, pr + se), 3)}
    agg = {}
    for s in systems:
        shu = [r for r in per if r["mode"] != "new" and r.get(s + "_auroc") is not None and r["n_pos"] >= 3]
        agg[s] = {"mean_auroc_shipped_units": round(float(np.mean([r[s + "_auroc"] for r in shu])), 3) if shu else None,
                  "calls_shipped_units": summ(pooled[s]["shipped_units"]),
                  "calls_new_units": summ(pooled[s]["new_units"]) if pooled[s]["new_units"].sum() else None,
                  "present_calls_per_patient_median": float(np.median([
                      sum(1 for r in systems[s][i].values() if r["call"] == "present") for i in range(len(rows))]))}
    return per, agg


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--expr", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--id-col", default="RhapsodyID")
    ap.add_argument("--exclude", default="WHO=Healthy control", help="COL=VALUE rows left out of training")
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", default=None)
    ap.add_argument("--sealed-frac", type=float, default=0.30)
    ap.add_argument("--nested-folds", type=int, default=5, help="0 skips the nested-CV estimate")
    ap.add_argument("--margin", type=float, default=0.03)
    ap.add_argument("--stability", type=float, default=0.70)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-jobs", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()

    z = np.load(a.expr, allow_pickle=True)
    X, genes, pts = z["X"].astype(float), [str(g) for g in z["genes"]], [str(p) for p in z["patients"]]
    tab = pd.read_csv(a.truth, sep="\t", low_memory=False).set_index(a.id_col).reindex(pts)
    col, val = a.exclude.split("=", 1) if a.exclude else (None, None)
    keep = ~tab[col].astype(str).eq(val).values if col else np.ones(len(pts), bool)
    excluded = np.where(~keep)[0]
    units = [dict(u) for u in UNITS_32266596]
    labels = {u["name"]: truth_vector(tab, u["truth"]) for u in units}
    for u in units:
        labels[u["name"]][~keep] = np.nan

    bp = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
    series = [pd.Series(X[i], index=genes) for i in range(len(pts))]
    print("shipped calls + decision scores for %d patients ..." % len(pts), flush=True)
    shipped = [bp.predict(s, ref="sc") for s in series]
    DEC = np.vstack([CA.shipped_decisions(bp, s, "sc") for s in series])

    rng = np.random.default_rng(a.seed)
    cand = np.where(keep)[0]
    sealed = np.sort(rng.choice(cand, int(round(a.sealed_frac * len(cand))), replace=False))
    dev = np.setdiff1d(cand, sealed)
    print("development %d | sealed %d | excluded %d" % (len(dev), len(sealed), len(excluded)), flush=True)

    def fit_on(rows, force_all=False):
        lab = {k: v[rows] for k, v in labels.items()}
        return CA.fit(bp, X[rows], genes, DEC[rows], lab, units, margin=a.margin, stability=a.stability, seed=a.seed,
                      n_jobs=a.n_jobs, shipped_calls=[shipped[i] for i in rows], force_all=force_all,
                      cohort=a.cohort)

    report = {"cohort": a.cohort, "n_development": int(len(dev)), "n_sealed": int(len(sealed)),
              "margin_auroc": a.margin}
    # ---- nested CV on the development set: the routing decision is refitted inside every outer fold ----
    if a.nested_folds:
        folds = np.random.default_rng(a.seed + 1).permutation(len(dev)) % a.nested_folds
        sysd = {"shipped": [None] * len(dev), "routed": [None] * len(dev), "cohort_model_everywhere": [None] * len(dev)}
        routed_counts = {}
        for k in range(a.nested_folds):
            tr, te = dev[folds != k], dev[folds == k]
            ad_r = fit_on(tr)
            ad_a = fit_on(tr, force_all=True)
            for u in ad_r.units:
                if u["routed"]:
                    routed_counts[u["name"]] = routed_counts.get(u["name"], 0) + 1
            for j, i in zip(np.where(folds == k)[0], te):
                sysd["shipped"][j] = shipped[i]
                sysd["routed"][j] = ad_r.apply(bp, series[i], shipped[i])
                sysd["cohort_model_everywhere"][j] = ad_a.apply(bp, series[i], shipped[i])
            print("  outer fold %d/%d done (%.0fs)" % (k + 1, a.nested_folds, time.time() - t0), flush=True)
        per, agg = score(units, labels, sysd, dev, fold=folds)
        report["nested_cv"] = {"per_unit": per, "summary": agg,
                               "times_routed_out_of_%d_folds" % a.nested_folds: routed_counts}
    # ---- final adapter on the whole development set, scored once on the sealed patients ----
    final = fit_on(dev)
    if len(sealed):
        sysd = {"shipped": [shipped[i] for i in sealed],
                "routed": [final.apply(bp, series[i], shipped[i]) for i in sealed]}
        per, agg = score(units, labels, sysd, sealed)
        report["sealed"] = {"per_unit": per, "summary": agg}
    if len(excluded):
        report["excluded_rows_present_calls"] = {
            "shipped": [sum(r["call"] == "present" for r in shipped[i].values()) for i in excluded],
            "routed": [sum(r["call"] == "present" for r in final.apply(bp, series[i], shipped[i]).values())
                       for i in excluded]}
    final.d["sealed_patients"] = [pts[i] for i in sealed]
    final.d["evaluation"] = {k: report[k]["summary"] for k in ("nested_cv", "sealed") if k in report}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    final.save(a.out)
    report["routing"] = final.summary()["routed"], final.summary()["kept_on_shipped"]
    if a.report:
        json.dump(report, open(a.report, "w"), indent=1, default=str)
    print(json.dumps({"routing": report["routing"], "nested_cv": report.get("nested_cv", {}).get("summary"),
                      "sealed": report.get("sealed", {}).get("summary")}, indent=1, default=str))
    print("wrote %s (%.0fs)" % (a.out, time.time() - t0))


if __name__ == "__main__":
    main()
