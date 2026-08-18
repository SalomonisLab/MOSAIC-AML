#!/usr/bin/env python3
"""Recalibrate CIPHER-AML's decision thresholds to a NEW cohort, using a small labelled subset.

This is the largest measured, unexploited gain in the platform, and it is a deployment protocol rather
than a modelling change. Measured by fitting per-category thresholds on half a cohort and transferring
them to the other half:

  cohort                 shipped F1   recalibrated F1   precision
  sealed single-cell        0.525          0.569        0.547 -> 0.596
  Trumpp/Waclawiczek        0.418          0.520        0.292 -> 0.542
  GSE281087                 0.250          0.407        0.182 -> 0.458
  Leucegene                 0.558          0.586        0.545 -> 0.651

Label-FREE recalibration does not work: within-cohort percentile matching and prevalence matching are
both WORSE than the shipped thresholds on three of those four cohorts. The gain requires ground truth
from the target cohort. Hence this tool: give it the specimens you have sequenced, and it returns
thresholds fitted for your cohort, with an honest estimate of what they buy.

WHAT IT WILL NOT DO. It refuses to emit a threshold for a category with too few labelled positives,
because a cut fitted on one or two positives is noise that will look like a calibration. Those
categories keep the shipped threshold and are listed as such.

  python calibrate_to_cohort.py --runs "runs/gse_*" --out labels/thresholds_gse281087.json
  python calibrate_to_cohort.py --runs "runs/gse_*" --holdout 0.5     # honest transfer estimate
"""
import os, sys, json, glob, argparse, collections, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MIN_POS = 3                      # labelled positives required before a category is recalibrated
MIN_SPECIMENS = 8


def load(pattern):
    """-> {category: [(probability, shipped_threshold, truth, specimen_index)]}"""
    by = collections.defaultdict(list)
    files = sorted(glob.glob(os.path.join(ROOT, pattern)))
    for si, f in enumerate(files):
        p = f if f.endswith(".json") else os.path.join(f, "patient_report.json")
        if not os.path.exists(p):
            continue
        try:
            rep = json.load(open(p, encoding="utf-8"))
        except Exception:
            continue
        for r in rep.get("mutation_predictions") or []:
            t = r.get("true_label")
            if t not in ("present", "absent") or r.get("probability") is None:
                continue
            by[r["mutation"]].append((float(r["probability"]), float(r.get("threshold") or 0.5),
                                      1 if t == "present" else 0, si))
    return by, len(files)


def counts(pred, y):
    return (int(np.sum((pred == 1) & (y == 1))), int(np.sum((pred == 1) & (y == 0))),
            int(np.sum((pred == 0) & (y == 1))), int(np.sum((pred == 0) & (y == 0))))


def f1_of(tp, fp, fn):
    pr = tp / (tp + fp) if tp + fp else 0.0
    sn = tp / (tp + fn) if tp + fn else 0.0
    return 2 * pr * sn / (pr + sn) if pr + sn > 0 else 0.0


def fit_thresholds(by, cats=None):
    """Per-category thresholds maximising the POOLED F1, by coordinate ascent.

    Pooled, not per-category: optimising F1 inside each category and then pooling is what the first
    version of this did, and it came out WORSE than the shipped thresholds on every cohort -- for a
    category with a single positive the F1-optimal cut is wildly permissive and those false positives
    pool. The metric optimised has to be the metric reported.
    """
    cats = cats or list(by)
    cur = {c: by[c][0][1] for c in cats}

    def pooled(cur):
        TP = FP = FN = TN = 0
        for c in cats:
            p = np.array([x[0] for x in by[c]]); y = np.array([x[2] for x in by[c]])
            a, b, d, e = counts((p >= cur[c]).astype(int), y)
            TP += a; FP += b; FN += d; TN += e
        return f1_of(TP, FP, FN), (TP, FP, FN, TN)

    best, _ = pooled(cur)
    for _ in range(6):
        improved = False
        for c in cats:
            p = np.array([x[0] for x in by[c]])
            keep = cur[c]
            for t in np.unique(np.concatenate([p, [1.01]])):
                cur[c] = float(t)
                f, _ = pooled(cur)
                if f > best + 1e-9:
                    best, keep, improved = f, float(t), True
            cur[c] = keep
        if not improved:
            break
    return cur, pooled(cur)


def summarise(by, thr):
    TP = FP = FN = TN = 0
    for c, v in by.items():
        p = np.array([x[0] for x in v]); y = np.array([x[2] for x in v])
        a, b, d, e = counts((p >= thr.get(c, v[0][1])).astype(int), y)
        TP += a; FP += b; FN += d; TN += e
    return {"tp": TP, "fp": FP, "fn": FN, "tn": TN,
            "sensitivity": round(TP / (TP + FN), 4) if TP + FN else None,
            "specificity": round(TN / (TN + FP), 4) if TN + FP else None,
            "precision": round(TP / (TP + FP), 4) if TP + FP else None,
            "f1": round(f1_of(TP, FP, FN), 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, help="glob of run directories or patient_report.json files")
    ap.add_argument("--out", default=None)
    ap.add_argument("--holdout", type=float, default=0.0,
                    help="fit on this fraction of specimens and score the rest, for an honest estimate")
    ap.add_argument("--name", default="cohort")
    a = ap.parse_args()

    by, n_files = load(a.runs)
    if not by:
        print("no labelled predictions matched %r" % a.runs); return
    specimens = sorted({x[3] for v in by.values() for x in v})
    print("%s: %d run(s), %d labelled specimens, %d categories, %d scored calls"
          % (a.name, n_files, len(specimens), len(by), sum(len(v) for v in by.values())))
    if len(specimens) < MIN_SPECIMENS:
        print("REFUSING: %d labelled specimens is below the minimum of %d; thresholds fitted on fewer "
              "would be noise." % (len(specimens), MIN_SPECIMENS))
        return

    eligible = [c for c, v in by.items() if sum(x[2] for x in v) >= MIN_POS]
    skipped = [c for c in by if c not in eligible]
    print("  %d categories have >= %d labelled positives and will be recalibrated; %d keep the shipped "
          "threshold" % (len(eligible), MIN_POS, len(skipped)))

    shipped = summarise(by, {c: by[c][0][1] for c in by})
    fitted, _ = fit_thresholds(by, eligible)
    full = {c: fitted.get(c, by[c][0][1]) for c in by}
    insample = summarise(by, full)

    out = {"cohort": a.name, "n_specimens": len(specimens), "n_categories": len(by),
           "recalibrated": len(eligible), "kept_shipped": sorted(skipped),
           "thresholds": {c: round(float(v), 4) for c, v in full.items()},
           "shipped_performance": shipped, "in_sample_performance": insample,
           "warning": ("in-sample figures are optimistic -- one threshold per category fitted on the "
                       "same labels it is scored against. Use --holdout for a transferable estimate.")}

    print("\n  %-22s %s" % ("shipped thresholds", json.dumps(shipped)))
    print("  %-22s %s" % ("recalibrated (in-sample)", json.dumps(insample)))

    if a.holdout > 0:
        rng = np.random.RandomState(0)
        perm = rng.permutation(specimens)
        k = int(round(a.holdout * len(specimens)))
        fit_set, score_set = set(perm[k:]), set(perm[:k])
        sub = {c: [x for x in v if x[3] in fit_set] for c, v in by.items()}
        sub = {c: v for c, v in sub.items() if v and sum(x[2] for x in v) >= MIN_POS}
        ft, _ = fit_thresholds(sub)
        held = {c: [x for x in v if x[3] in score_set] for c, v in by.items()}
        held = {c: v for c, v in held.items() if v}
        tr = summarise(held, {c: ft.get(c, by[c][0][1]) for c in held})
        sh = summarise(held, {c: by[c][0][1] for c in held})
        out["holdout_estimate"] = {"fit_specimens": len(fit_set), "scored_specimens": len(score_set),
                                   "shipped": sh, "recalibrated": tr,
                                   "delta_f1": round((tr["f1"] or 0) - (sh["f1"] or 0), 4)}
        print("\n  transfer estimate (fit on %d specimens, scored on %d):"
              % (len(fit_set), len(score_set)))
        print("    shipped      F1 %s  precision %s" % (sh["f1"], sh["precision"]))
        print("    recalibrated F1 %s  precision %s   (delta %+.4f)"
              % (tr["f1"], tr["precision"], out["holdout_estimate"]["delta_f1"]))

    dst = a.out or os.path.join(ROOT, "labels", "thresholds_%s.json" % a.name)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    json.dump(out, open(dst, "w"), indent=1)
    print("\n  wrote %s" % dst)


if __name__ == "__main__":
    main()
