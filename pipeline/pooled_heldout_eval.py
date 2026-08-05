#!/usr/bin/env python3
"""POOLED held-out evaluation — every mutation counted, no per-category minimum.

NS's critique was correct: the earlier "held-out scRNA" row reported only 3 mutations because that
analysis required >=3 positives AND >=3 negatives per category before scoring it. For an OVERALL
sensitivity/specificity that filter is wrong — a mutation seen once still contributes one true
positive or one false negative. This pools every (sample x mutation) call with a known label across
ALL held-out single-cell cohorts and micro-averages.

  cohorts: sealed held-out scRNA (runs/predict_*, n=29) + Trumpp/Waclawiczek (runs/trumpp_*, n=16)
  reported: pooled TP/FP/FN/TN -> overall sensitivity, specificity, precision, F1
            plus per-mutation rows for ANY mutation with >=1 positive (no minimum)

  /usr/local/anaconda3-2020/bin/python pooled_heldout_eval.py     (runs locally too)
"""
import os, json, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)

def canon(c):
    cl = str(c).lower()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl: return "inv16"
    if "kmt2a" in cl: return "KMT2A"
    return str(c)

def collect(pattern, label):
    """-> {mutation: [tp,fp,fn,tn]}, n_samples, n_calls"""
    agg, ns, nc = {}, 0, 0
    for f in sorted(glob.glob(os.path.join(ROOT, "runs", pattern, "patient_report.json"))):
        rep = json.load(open(f)); ns += 1
        preds = rep.get("mutation_predictions") or []
        known = set(str(x).upper() for x in (rep.get("known_drivers") or []))
        has_known = bool(rep.get("known_drivers"))
        for p in preds:
            m = canon(p["mutation"]); call = (p.get("call") == "present")
            tl = p.get("true_label")
            if tl is not None:
                truth = (tl == "present")
            elif has_known:
                # Table S4-style cohorts list the drivers; anything else the model scores is wild-type
                truth = str(p["mutation"]).split("_")[0].split("-")[0].upper() in known or str(m).upper() in known
            else:
                continue
            d = agg.setdefault(m, [0, 0, 0, 0]); nc += 1
            d[0 if (truth and call) else 2 if truth else 1 if call else 3] += 1
    return agg, ns, nc

def metrics(tp, fp, fn, tn):
    se = tp / (tp + fn) if tp + fn else None
    sp = tn / (tn + fp) if tn + fp else None
    pr = tp / (tp + fp) if tp + fp else None
    f1 = (2 * pr * se / (pr + se)) if (pr and se) else None
    acc = (tp + tn) / max(tp + fp + fn + tn, 1)
    return dict(tp=tp, fp=fp, fn=fn, tn=tn,
                sensitivity=round(se, 4) if se is not None else None,
                specificity=round(sp, 4) if sp is not None else None,
                precision=round(pr, 4) if pr is not None else None,
                f1=round(f1, 4) if f1 is not None else None, accuracy=round(acc, 4))

COH = [("predict_*", "sealed held-out scRNA"), ("trumpp_*", "Trumpp/Waclawiczek"),
       ("gse_*", "GSE281087 (panel-honest)")]
out, grand = {}, [0, 0, 0, 0]
for pat, lab in COH:
    agg, ns, nc = collect(pat, lab)
    tot = [sum(v[i] for v in agg.values()) for i in range(4)]
    for i in range(4): grand[i] += tot[i]
    permut = {m: metrics(*v) for m, v in sorted(agg.items(), key=lambda kv: -(kv[1][0] + kv[1][2]))}
    npos_any = sum(1 for v in agg.values() if v[0] + v[2] >= 1)
    out[lab] = {"n_samples": ns, "n_calls_scored": nc, "n_mutations_any_positive": npos_any,
                "n_mutations_ge3_positive": sum(1 for v in agg.values() if v[0] + v[2] >= 3),
                "overall": metrics(*tot), "per_mutation": permut}
    o = out[lab]["overall"]
    print("=== %s: %d samples, %d calls scored ===" % (lab, ns, nc))
    print("    mutations with >=1 positive: %d   (with >=3: %d)" % (npos_any, out[lab]["n_mutations_ge3_positive"]))
    print("    POOLED sensitivity %s  specificity %s  precision %s  F1 %s  (TP=%d FP=%d FN=%d TN=%d)"
          % (o["sensitivity"], o["specificity"], o["precision"], o["f1"], o["tp"], o["fp"], o["fn"], o["tn"]))

g = metrics(*grand)
out["ALL held-out single-cell"] = {"overall": g,
    "n_samples": sum(out[l]["n_samples"] for _, l in COH),
    "n_calls_scored": sum(out[l]["n_calls_scored"] for _, l in COH)}
print("\n=== ALL held-out single-cell pooled (%d samples, %d calls) ===" % (
    out["ALL held-out single-cell"]["n_samples"], out["ALL held-out single-cell"]["n_calls_scored"]))
print("    sensitivity %s   specificity %s   precision %s   F1 %s   accuracy %s"
      % (g["sensitivity"], g["specificity"], g["precision"], g["f1"], g["accuracy"]))
print("    TP=%d FP=%d FN=%d TN=%d" % (g["tp"], g["fp"], g["fn"], g["tn"]))

# per-mutation table (>=1 positive, no minimum) for the sealed set
print("\n=== per-mutation, sealed held-out (>=1 positive, NO minimum filter) ===")
pm = out["sealed held-out scRNA"]["per_mutation"]
print("%-22s %4s %4s %4s %4s  %6s %6s" % ("mutation", "TP", "FP", "FN", "TN", "sens", "spec"))
for m, v in pm.items():
    if v["tp"] + v["fn"] >= 1:
        print("%-22s %4d %4d %4d %4d  %6s %6s" % (m, v["tp"], v["fp"], v["fn"], v["tn"], v["sensitivity"], v["specificity"]))

json.dump(out, open(os.path.join(OUTD, "pooled_heldout_eval.json"), "w"), indent=1)
print("\nwrote deliverables/pooled_heldout_eval.json")
print("POOLED HELDOUT OK")
