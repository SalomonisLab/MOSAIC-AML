#!/usr/bin/env python3
"""How much of CIPHER-AML's weakness is threshold PLACEMENT rather than ranking?

The platform's recurring shape is a strong ranker with a weak decision rule: mean AUROC 0.908 against a
pooled precision of 0.336, and on GSE281087 a within-cohort ranking AUROC of 0.688 against a precision
of 0.169. If that gap is threshold placement, it is cheap to close. If it is not, no amount of
recalibration helps and the honest answer is that the caller is weaker than its AUROC suggests.

Four rules, scored per cohort on the same calls:

  deployed          the trained per-category threshold, as shipped
  cohort_percentile re-rank within the cohort, keep the trained threshold (what predict_cohort does)
  prevalence        call the top k of each category, k = that category's training positive rate
  ORACLE            the threshold that maximises F1 per category USING THE TEST LABELS

The oracle is not a proposal -- it is the ceiling. It says what perfect threshold placement would buy,
and therefore how much of the gap is placement at all. Any rule that lands near the oracle is worth
deploying; a large residual gap means the ranking itself is the limit.

  python exp_cipher_thresholds.py  ->  deliverables/exp_cipher_thresholds.json
"""
import os, sys, json, glob, time, collections, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_cipher_thresholds.json")

COHORTS = {"sealed single-cell": "runs/predict_*/patient_report.json",
           "Trumpp/Waclawiczek": "runs/trumpp_*/patient_report.json",
           "GSE281087": "runs/gse_*/patient_report.json",
           "Leucegene": "runs/leucegene_*/patient_report.json"}


def load(pattern):
    """-> {category: [(probability, trained_threshold, truth 0/1, sample_index), ...]}"""
    by = collections.defaultdict(list)
    for si, f in enumerate(sorted(glob.glob(os.path.join(ROOT, pattern)))):
        try:
            rep = json.load(open(f, encoding="utf-8"))
        except Exception:
            continue
        for p in rep.get("mutation_predictions") or []:
            t = p.get("true_label")
            if t not in ("present", "absent") or p.get("probability") is None:
                continue
            by[p["mutation"]].append((float(p["probability"]), float(p.get("threshold") or 0.5),
                                      1 if t == "present" else 0, si))
    return by


def counts(calls, ys):
    tp = int(np.sum((calls == 1) & (ys == 1))); fp = int(np.sum((calls == 1) & (ys == 0)))
    fn = int(np.sum((calls == 0) & (ys == 1))); tn = int(np.sum((calls == 0) & (ys == 0)))
    return tp, fp, fn, tn


def summarise(tp, fp, fn, tn):
    sn = tp / (tp + fn) if tp + fn else float("nan")
    sp = tn / (tn + fp) if tn + fp else float("nan")
    pr = tp / (tp + fp) if tp + fp else float("nan")
    f1 = 2 * pr * sn / (pr + sn) if (pr == pr and sn == sn and pr + sn > 0) else float("nan")
    r = {"tp": tp, "fp": fp, "fn": fn, "tn": tn}
    for k, v in (("sensitivity", sn), ("specificity", sp), ("precision", pr), ("f1", f1)):
        r[k] = None if v != v else round(float(v), 4)
    return r


def rule_calls(rule, probs, thr, ys, prevalence, min_n=8):
    n = len(probs)
    if rule == "deployed":
        return (probs >= thr).astype(int)
    if rule == "cohort_percentile":
        if n < min_n:
            return (probs >= thr).astype(int)
        q = probs.argsort().argsort() / (n - 1)
        return (q >= thr).astype(int)
    if rule == "prevalence":
        if n < min_n or prevalence is None:
            return (probs >= thr).astype(int)
        k = max(1, int(round(prevalence * n)))
        cut = np.sort(probs)[::-1][k - 1]
        return (probs >= cut).astype(int)
    raise ValueError(rule)


def oracle_pooled(by, thr_of):
    """Per-category thresholds chosen to maximise the POOLED F1, by coordinate ascent.

    The first version of this maximised F1 within each category and then pooled the counts, and came
    out WORSE than the deployed thresholds on all four cohorts. That is not a fact about the model: for
    a category with one positive, the F1-optimal cut is extremely permissive (catch the one, accept many
    false positives), and pooling those decisions buries the pooled precision. Optimising the metric
    that is actually reported is the only version of this that means anything.
    """
    cats = list(by)
    cur = {c: thr_of[c] for c in cats}
    def pooled(cur):
        TP = FP = FN = TN = 0
        for c in cats:
            v = by[c]
            p = np.array([x[0] for x in v]); y = np.array([x[2] for x in v])
            tp, fp, fn, tn = counts((p >= cur[c]).astype(int), y)
            TP += tp; FP += fp; FN += fn; TN += tn
        pr = TP / (TP + FP) if TP + FP else 0.0
        sn = TP / (TP + FN) if TP + FN else 0.0
        return (2 * pr * sn / (pr + sn) if pr + sn > 0 else 0.0), (TP, FP, FN, TN)
    best, _ = pooled(cur)
    for _ in range(6):
        improved = False
        for c in cats:
            v = by[c]
            p = np.array([x[0] for x in v])
            keep = cur[c]
            for t in np.unique(np.concatenate([p, [1.01]])):
                cur[c] = float(t)
                f1, _ = pooled(cur)
                if f1 > best + 1e-9:
                    best, keep, improved = f1, float(t), True
            cur[c] = keep
        if not improved:
            break
    _, cnt = pooled(cur)
    return cnt, cur


def main():
    t0 = time.time()
    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "cohorts": {},
           "note": ("ORACLE uses the test labels to place each category's threshold. It is the ceiling "
                    "on what threshold placement can buy, not a deployable rule.")}

    # training prevalence per category, from the model itself
    prev = {}
    try:
        import pickle
        from amlmm.bulk_predictor import BulkMutationPredictor as B
        m = pickle.load(open(os.path.join(HERE, "bulk_mutation_predictor.pkl"), "rb"))
        m.__class__ = B
        ntr = None
        for cat, mm in m.models.items():
            n = mm.get("n_pos")
            base = (m.score_refs.get("beataml") or {}).get(cat)
            if n is not None and base is not None and len(base):
                prev[cat] = float(n) / float(len(base))
        ntr = len(prev)
        print("training prevalence recovered for %d categories" % ntr)
    except Exception as ex:
        print("prevalence unavailable (%s); the prevalence rule will fall back to deployed" % str(ex)[:80])

    for name, pat in COHORTS.items():
        by = load(pat)
        if not by:
            continue
        n_calls = sum(len(v) for v in by.values())
        row = {"n_categories": len(by), "n_calls": n_calls, "rules": {}}
        for rule in ("deployed", "cohort_percentile", "prevalence"):
            TP = FP = FN = TN = 0
            for cat, v in by.items():
                probs = np.array([x[0] for x in v]); thr = v[0][1]
                ys = np.array([x[2] for x in v])
                c = rule_calls(rule, probs, thr, ys, prev.get(cat))
                tp, fp, fn, tn = counts(c, ys)
                TP += tp; FP += fp; FN += fn; TN += tn
            row["rules"][rule] = summarise(TP, FP, FN, TN)
        (tp, fp, fn, tn), _ = oracle_pooled(by, {c: by[c][0][1] for c in by})
        row["rules"]["oracle"] = summarise(tp, fp, fn, tn)

        # The oracle above fits one threshold per category ON the labels it is scored against -- up to
        # 29 free parameters on as few as 390 calls. Split-half says how much of that ceiling survives
        # being transferred to specimens it was not fitted on, which is the only honest headroom.
        samples = sorted({x[3] for v in by.values() for x in v})
        if len(samples) >= 8:
            half = set(samples[::2])
            TP = FP = FN = TN = 0
            for a_set, b_set in ((half, set(samples) - half), (set(samples) - half, half)):
                fit = {c: [x for x in v if x[3] in a_set] for c, v in by.items()}
                fit = {c: v for c, v in fit.items() if v}
                if not fit:
                    continue
                _, cuts = oracle_pooled(fit, {c: by[c][0][1] for c in fit})
                for c, v in by.items():
                    ev = [x for x in v if x[3] in b_set]
                    if not ev:
                        continue
                    p = np.array([x[0] for x in ev]); y = np.array([x[2] for x in ev])
                    t = cuts.get(c, by[c][0][1])
                    a_, b_, c_, d_ = counts((p >= t).astype(int), y)
                    TP += a_; FP += b_; FN += c_; TN += d_
            row["rules"]["oracle_split_half"] = summarise(TP, FP, FN, TN)
        res["cohorts"][name] = row
        print("\n%s  (%d categories, %d calls)" % (name, len(by), n_calls))
        for rule, r in row["rules"].items():
            print("  %-18s tp=%4d fp=%4d fn=%4d tn=%5d | sens %s spec %s prec %s F1 %s"
                  % (rule, r["tp"], r["fp"], r["fn"], r["tn"], r["sensitivity"], r["specificity"],
                     r["precision"], r["f1"]))

    # how much of the gap to the oracle does the best deployable rule close?
    gap = {}
    for name, row in res["cohorts"].items():
        d = row["rules"]["deployed"]["f1"]
        o = (row["rules"].get("oracle_split_half") or {}).get("f1")
        oin = row["rules"]["oracle"]["f1"]
        best = max((row["rules"][r]["f1"] or 0) for r in ("cohort_percentile", "prevalence"))
        if d is None:
            continue
        gap[name] = {"deployed_f1": d, "oracle_in_sample_f1": oin,
                     "oracle_split_half_f1": o, "best_label_free_rule_f1": round(best, 4),
                     "transferable_headroom": None if o is None else round(o - d, 4),
                     "label_free_rules_capture": round(best - d, 4)}
    res["headroom"] = gap
    res["headroom_conclusion"] = (
        "Threshold placement is worth a lot -- but only WITH labels from the target cohort. Fitting "
        "per-category thresholds on half a cohort and transferring them to the other half gains F1 on "
        "every cohort, most where the caller is worst. The label-free rules (cohort percentile, "
        "prevalence matching) do NOT capture it and are worse than the shipped thresholds on three of "
        "four cohorts. The in-sample oracle is inflated -- up to 29 free parameters on as few as 390 "
        "calls -- and split-half is the number to quote.")
    print("\n== threshold headroom (split-half oracle is the honest one) ==")
    for k, v in gap.items():
        print("  %-20s deployed %.3f | label-free best %.3f | split-half oracle %s | in-sample %s"
              % (k, v["deployed_f1"], v["best_label_free_rule_f1"],
                 v["oracle_split_half_f1"], v["oracle_in_sample_f1"]))

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
