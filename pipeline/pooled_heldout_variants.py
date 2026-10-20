#!/usr/bin/env python3
"""How the pooled held-out figure depends on where its labels come from.

pooled_heldout_eval.py scores every (sample x mutation) call with a label. Where a cohort
lists its drivers but records no label for a call, it infers one: positive if the call's gene
prefix is a listed driver, negative otherwise. This script holds everything else fixed and
varies one choice at a time, so each effect can be read on its own:

  labels   L0  as shipped: recorded, and inferred where none is recorded
           L1  recorded labels only
           L2  recorded, plus inferred negatives (inferred positives dropped)
           L3  labels at the resolution the truth has: recorded labels; an unrecorded
               subtype takes the label of a recorded call naming the same lesion; a listed
               driver whose subtype is not recorded becomes one gene-level positive,
               "did the model call any subtype of it", rather than a positive per subtype
  counting D   one call per patient and lesion (aliases and second model heads once)
           C   controls excluded

It also tests the inference rule itself on the calls where the truth is recorded.

Patient-level bootstrap (2000 resamples, seed 20261020) gives 95% intervals, and paired
differences between variants use the same resamples. Results are reported per model:
the sealed held-out calls come from the multimodal predictor, the Trumpp/Waclawiczek and
GSE281087 calls from the bulk caller run on single-cell input.

Writes deliverables/pooled_heldout_variants.json and .tsv. Runs locally, standard library only.
"""
import glob
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import categories as cat

OUTD = os.path.join(ROOT, "deliverables")
SEED, B = 20261020, 2000
COHORTS = [("predict_*", "sealed held-out scRNA", "multimodal"),
           ("trumpp_*", "Trumpp/Waclawiczek", "bulk caller"),
           ("gse_*", "GSE281087 (panel-honest)", "bulk caller")]


def shipped_rule(name, known):
    """pooled_heldout_eval.py's inferred label, reproduced exactly."""
    cl = str(name).lower()
    canon = "inv16" if ("inv16" in cl or "inv(16)" in cl or "cbfb" in cl) else \
        "KMT2A" if "kmt2a" in cl else str(name)
    return str(name).split("_")[0].split("-")[0].upper() in known or canon.upper() in known


def shipped_key(name):
    cl = str(name).lower()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl:
        return "inv16"
    if "kmt2a" in cl:
        return "KMT2A"
    return str(name)


def lesion(name):
    """One key per lesion: aliases and a second model head share it; a bulk
    category that names a multimodal category's lesion shares it too."""
    m, how = cat.match(name)
    if m and how == "same":
        name = m
    name = cat.row(name)
    return cat.SECOND_HEADS.get(name, name)


def load():
    patients = []
    for pat, cohort, model in COHORTS:
        for f in sorted(glob.glob(os.path.join(ROOT, "runs", pat, "patient_report.json"))):
            r = json.load(open(f, encoding="utf-8"))
            known = set(str(x).upper() for x in (r.get("known_drivers") or []))
            calls = [{"name": p["mutation"], "call": p.get("call") == "present",
                      "tl": p.get("true_label")} for p in (r.get("mutation_predictions") or [])]
            patients.append({"id": os.path.basename(os.path.dirname(f)), "cohort": cohort,
                             "model": model, "known": known, "has_known": bool(r.get("known_drivers")),
                             "control": r.get("specimen_class") == "control", "calls": calls})
    return patients


def units(p, labels, dedupe):
    """(key, truth, call) units one patient contributes under a labeling policy."""
    out = []
    recorded = {c["name"]: c["tl"] == "present" for c in p["calls"] if c["tl"] is not None}
    rec_lesion = {}
    for n, t in recorded.items():
        rec_lesion.setdefault(lesion(n), t)
    gene_level = {}
    for c in p["calls"]:
        n, call = c["name"], c["call"]
        if c["tl"] is not None:
            out.append((n, c["tl"] == "present", call))
            continue
        if not p["has_known"] or labels == "L1":
            continue
        if labels in ("L0", "L2"):
            truth = shipped_rule(n, p["known"])
            if labels == "L2" and truth:
                continue
            out.append((n, truth, call))
            continue
        # L3
        if lesion(n) in rec_lesion:
            out.append((n, rec_lesion[lesion(n)], call))
            continue
        g = cat.gene(n)
        if g.upper() not in p["known"]:
            out.append((n, False, call))
            continue
        if any(t for k, t in recorded.items() if cat.gene(k) == g):
            continue            # a recorded subtype explains the driver; this one is unknown
        gene_level[g] = gene_level.get(g, False) or call
    for g, call in gene_level.items():
        out.append(("%s (gene level)" % g, True, call))
    if dedupe:
        seen, kept = {}, []
        for k, t, call in out:
            key = lesion(k) if dedupe == "lesion" else shipped_key(k)
            if key in seen:
                continue
            seen[key] = True
            kept.append((k, t, call))
        out = kept
    return out


def summarize(pts, labels, dedupe, controls):
    tp = fp = fn = tn = 0
    per_cat, per_pat = {}, []
    for p in pts:
        if controls == "excluded" and p["control"]:
            continue
        u = units(p, labels, dedupe)
        ppos = pcall = 0
        for k, t, c in u:
            if t and c:
                tp += 1
            elif t:
                fn += 1
            elif c:
                fp += 1
            else:
                tn += 1
            d = per_cat.setdefault(k, [0, 0, 0, 0])
            d[0 if (t and c) else 2 if t else 1 if c else 3] += 1
            if t:
                ppos += 1
                pcall += c
        if ppos:
            per_pat.append(pcall / ppos)
    se = tp / (tp + fn) if tp + fn else None
    sp = tn / (tn + fp) if tn + fp else None
    pr = tp / (tp + fp) if tp + fp else None
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else None
    cats = [v for v in per_cat.values() if v[0] + v[2]]
    macro = sum(v[0] / (v[0] + v[2]) for v in cats) / len(cats) if cats else None
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "sensitivity": se, "specificity": sp,
            "precision": pr, "f1": f1, "macro_sensitivity_per_category": macro,
            "n_categories_with_positives": len(cats),
            "mean_sensitivity_per_patient": (sum(per_pat) / len(per_pat)) if per_pat else None,
            "n_patients_with_positives": len(per_pat)}


def rule_check(pts):
    """The inference rule, tested on the calls whose truth is recorded."""
    out = {"TP": 0, "FP": 0, "FN": 0, "TN": 0, "disagreements": {}}
    for p in pts:
        if not p["has_known"]:
            continue
        for c in p["calls"]:
            if c["tl"] is None:
                continue
            truth, guess = c["tl"] == "present", shipped_rule(c["name"], p["known"])
            key = "TP" if truth and guess else "FN" if truth else "FP" if guess else "TN"
            out[key] += 1
            if key in ("FP", "FN"):
                out["disagreements"].setdefault(c["name"], {"FP": 0, "FN": 0})[key] += 1
    n_pos_rule = out["TP"] + out["FP"]
    out["share_of_rule_positives_wrong"] = out["FP"] / n_pos_rule if n_pos_rule else None
    n_neg_rule = out["TN"] + out["FN"]
    out["share_of_rule_negatives_wrong"] = out["FN"] / n_neg_rule if n_neg_rule else None
    return out


VARIANTS = [("L0", None, "included"), ("L1", None, "included"), ("L2", None, "included"),
            ("L3", None, "included"), ("L0", "lesion", "excluded"), ("L1", "lesion", "excluded"),
            ("L3", "lesion", "excluded")]
NAMES = {("L0", None, "included"): "as shipped",
         ("L1", None, "included"): "recorded labels only",
         ("L2", None, "included"): "recorded + inferred negatives",
         ("L3", None, "included"): "resolution-matched",
         ("L0", "lesion", "excluded"): "as shipped, one call per lesion, no controls",
         ("L1", "lesion", "excluded"): "recorded only, one call per lesion, no controls",
         ("L3", "lesion", "excluded"): "resolution-matched, one call per lesion, no controls"}
METRICS = ("sensitivity", "specificity", "precision", "f1", "macro_sensitivity_per_category",
           "mean_sensitivity_per_patient")


def pct(xs, q):
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    i = q * (len(xs) - 1)
    lo, hi = int(i), min(int(i) + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)


def analyze(pts):
    point = {v: summarize(pts, *v) for v in VARIANTS}
    rng = random.Random(SEED)
    boots = {v: {m: [] for m in METRICS} for v in VARIANTS}
    for _ in range(B):
        sample = [pts[rng.randrange(len(pts))] for _ in pts]
        for v in VARIANTS:
            s = summarize(sample, *v)
            for m in METRICS:
                boots[v][m].append(s[m])
    res = {}
    for v in VARIANTS:
        r = dict(point[v])
        for m in METRICS:
            r[m + "_ci95"] = [pct(boots[v][m], 0.025), pct(boots[v][m], 0.975)]
        res[NAMES[v]] = r
    # paired differences against the shipped labels, on the same resamples
    base = VARIANTS[0]
    diffs = {}
    for v in VARIANTS[1:]:
        d = {}
        for m in ("sensitivity", "specificity"):
            ds = [a - b for a, b in zip(boots[v][m], boots[base][m]) if a is not None and b is not None]
            pe = (point[v][m] - point[base][m]) if point[v][m] is not None and point[base][m] is not None else None
            d[m] = {"difference": pe, "ci95": [pct(ds, 0.025), pct(ds, 0.975)]}
        diffs[NAMES[v] + " - as shipped"] = d
    return res, diffs


def main():
    pts = load()
    out = {"seed": SEED, "bootstrap_resamples": B, "cohorts": {}, "models": {}}
    shipped = json.load(open(os.path.join(OUTD, "pooled_heldout_eval.json"), encoding="utf-8"))
    for _pat, cohort, model in COHORTS:
        cp = [p for p in pts if p["cohort"] == cohort]
        res, diffs = analyze(cp)
        s = res["as shipped"]
        ship = shipped[cohort]["overall"]
        if (s["tp"], s["fp"], s["fn"], s["tn"]) != (ship["tp"], ship["fp"], ship["fn"], ship["tn"]):
            raise SystemExit("%s: the shipped labels do not reproduce pooled_heldout_eval.json" % cohort)
        out["cohorts"][cohort] = {"model": model, "n_patients": len(cp), "variants": res,
                                  "paired_differences": diffs, "rule_check": rule_check(cp)}
    for model in ("multimodal", "bulk caller"):
        mp = [p for p in pts if p["model"] == model]
        res, diffs = analyze(mp)
        out["models"][model] = {"cohorts": sorted({p["cohort"] for p in mp}), "n_patients": len(mp),
                                "variants": res, "paired_differences": diffs, "rule_check": rule_check(mp)}
    res, diffs = analyze(pts)
    out["all pooled, both models"] = {"n_patients": len(pts), "variants": res, "paired_differences": diffs}
    json.dump(out, open(os.path.join(OUTD, "pooled_heldout_variants.json"), "w", encoding="utf-8"), indent=1)

    def f(x):
        return "" if x is None else "%.3f" % x
    with open(os.path.join(OUTD, "pooled_heldout_variants.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("scope\tmodel\tvariant\tn_patients\ttp\tfp\tfn\ttn\tsensitivity\tsens_lo\tsens_hi\t"
                 "specificity\tspec_lo\tspec_hi\tprecision\tf1\tmacro_sensitivity_per_category\t"
                 "mean_sensitivity_per_patient\n")
        blocks = [(c, v["model"], v) for c, v in out["cohorts"].items()] + \
                 [("model: " + m, m, v) for m, v in out["models"].items()] + \
                 [("all pooled, both models", "both", out["all pooled, both models"])]
        for scope, model, blk in blocks:
            for name, r in blk["variants"].items():
                fh.write("\t".join([scope, model, name, str(blk["n_patients"]), str(r["tp"]), str(r["fp"]),
                                    str(r["fn"]), str(r["tn"]), f(r["sensitivity"]), f(r["sensitivity_ci95"][0]),
                                    f(r["sensitivity_ci95"][1]), f(r["specificity"]), f(r["specificity_ci95"][0]),
                                    f(r["specificity_ci95"][1]), f(r["precision"]), f(r["f1"]),
                                    f(r["macro_sensitivity_per_category"]),
                                    f(r["mean_sensitivity_per_patient"])]) + "\n")
    for scope, blk in list(out["cohorts"].items()) + [("model: " + m, v) for m, v in out["models"].items()]:
        print("=== %s (%d patients)" % (scope, blk["n_patients"]))
        for name, r in blk["variants"].items():
            print("  %-52s sens %s [%s, %s]  spec %s  TP/FP/FN/TN %d/%d/%d/%d" % (
                name, f(r["sensitivity"]), f(r["sensitivity_ci95"][0]), f(r["sensitivity_ci95"][1]),
                f(r["specificity"]), r["tp"], r["fp"], r["fn"], r["tn"]))
        rc = blk["rule_check"]
        if rc["TP"] + rc["FP"] + rc["FN"] + rc["TN"]:
            print("  inference rule on recorded calls: TP %d FP %d FN %d TN %d" % (rc["TP"], rc["FP"], rc["FN"], rc["TN"]))
    print("wrote deliverables/pooled_heldout_variants.json and .tsv")
    print("POOLED VARIANTS OK")


if __name__ == "__main__":
    main()
