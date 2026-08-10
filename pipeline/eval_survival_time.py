#!/usr/bin/env python3
"""The second half of the question: not "who dies first" but "how long, and can we believe the number".

C-index answers only the ranking question. A model can rank patients perfectly and still be badly wrong
about *when*. This script asks the harder thing three ways, all on the sealed hold-out:

  1. **Calibration at a horizon.** Bin patients by predicted P(alive at 1/2/5 y), compare the mean
     prediction in each bin with the Kaplan-Meier observed survival of that bin. Predicted 70% should
     mean 70% alive.
  2. **Point-estimate error on median survival.** For each risk tertile, the model's median predicted
     survival against the group's observed Kaplan-Meier median. A per-patient median is also produced,
     but the error on it is what decides whether it is honest to show one.
  3. **How wide the honest interval is.** The spread of actual survival times inside a predicted-risk
     band — the real answer to "how long will this patient live", and the reason a single number is
     misleading even when the ranking is good.

  python eval_survival_time.py -> deliverables/survival_time_validation.json
"""
import os, sys, json, pickle, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold

from amlmm.survival import data as SD, coxph as CX
from train_survival_model import build_blocks, fit_arm, ARM_BLOCKS, HORIZONS, BUNDLE

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "survival_time_validation.json")


def km_at(t, e, horizon):
    ts, ss = CX.km(t, e)
    if not len(ts):
        return 1.0
    return float(np.interp(horizon, ts, ss, left=1.0, right=ss[-1]))


def km_median(t, e):
    ts, ss = CX.km(t, e)
    below = np.where(ss <= 0.5)[0]
    return float(ts[below[0]]) if len(below) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--holdout", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-pc", type=int, default=60)
    a = ap.parse_args()

    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(x) for x in d["genes"]]
    ba = [str(s) for s in d["ba_samples"]]
    X_lin = d["ba_X"].astype(np.float64)
    mut_all = pd.DataFrame(d["ba_L"].astype(float), index=ba, columns=[str(c) for c in d["drivers"]])
    row_of = {s: i for i, s in enumerate(ba)}

    cl = SD.load_cohort(specimens=ba)
    rows = np.array([row_of[s] for s in cl["specimen"]])
    t = cl["time_years"].values.astype(float)
    e = cl["event"].values.astype(int)
    g = cl["subject"].values

    subs = np.array(sorted(set(g)))
    rng = np.random.RandomState(a.seed)
    perm = rng.permutation(len(subs))
    hold = set(subs[perm[:int(round(a.holdout * len(subs)))]])
    is_hold = np.array([s in hold for s in g]); tr = ~is_hold

    fs, B = build_blocks(cl, X_lin, mut_all, rows, sorted(rows[tr]), genes, a.n_pc, 4000)
    Btr = {k: v[tr] for k, v in B.items()}
    Bho = {k: v[is_hold] for k, v in B.items()}
    r_ho, mod = fit_arm("full", Btr, Bho, t[tr], e[tr], g[tr])
    S_ho = mod.survival({b: Bho[b] for b in ARM_BLOCKS["full"]}, HORIZONS)
    med_ho = mod.median_survival({b: Bho[b] for b in ARM_BLOCKS["full"]})
    th, eh = t[is_hold], e[is_hold]

    res = {"n_holdout": int(is_hold.sum()), "events_holdout": int(eh.sum()), "horizons": HORIZONS}

    # ---- 1. calibration at each horizon ----
    cal = {}
    for j, h in enumerate(HORIZONS):
        p = S_ho[:, j]
        q = np.quantile(p, [0, .25, .5, .75, 1.0])
        bins = []
        for i in range(4):
            m = (p >= q[i]) & (p <= q[i + 1] if i == 3 else p < q[i + 1])
            if m.sum() < 8:
                continue
            bins.append({"n": int(m.sum()), "predicted": round(float(p[m].mean()), 3),
                         "observed_km": round(km_at(th[m], eh[m], h), 3)})
        gap = [abs(b["predicted"] - b["observed_km"]) for b in bins]
        cal["%gy" % h] = {"bins": bins,
                          "mean_abs_calibration_gap": round(float(np.mean(gap)), 3) if gap else None,
                          "overall_predicted": round(float(S_ho[:, j].mean()), 3),
                          "overall_observed_km": round(km_at(th, eh, h), 3)}
    res["calibration"] = cal

    # ---- 2. median survival, per risk tertile ----
    cut = np.quantile(r_ho, [1 / 3, 2 / 3])
    grp = np.digitize(r_ho, cut)
    tert = {}
    for gi, name in enumerate(("low", "intermediate", "high")):
        m = grp == gi
        if m.sum() < 5:
            continue
        obs = km_median(th[m], eh[m])
        pred = float(np.nanmedian(med_ho[m])) if np.isfinite(med_ho[m]).any() else None
        tert[name] = {
            "n": int(m.sum()), "events": int(eh[m].sum()),
            "predicted_median_years": None if pred is None else round(pred, 2),
            "observed_km_median_years": None if obs is None else round(obs, 2),
            "error_years": None if (pred is None or obs is None) else round(abs(pred - obs), 2),
            "observed_survival_time_iqr_years": [round(float(np.percentile(th[m & (eh == 1)], 25)), 2),
                                                 round(float(np.percentile(th[m & (eh == 1)], 75)), 2)]
                                                if (m & (eh == 1)).sum() >= 4 else None,
        }
    res["median_survival_by_risk_tertile"] = tert

    # ---- 3. how wide is the honest interval for one patient ----
    died = eh == 1
    spread = {}
    for gi, name in enumerate(("low", "intermediate", "high")):
        m = (grp == gi) & died
        if m.sum() >= 6:
            spread[name] = {"n_deaths": int(m.sum()),
                            "survival_years_p10_p90": [round(float(np.percentile(th[m], 10)), 2),
                                                       round(float(np.percentile(th[m], 90)), 2)],
                            "median": round(float(np.median(th[m])), 2)}
    res["individual_spread_within_risk_group"] = spread

    # a blunt summary of how far a per-patient point estimate is from the truth, among those who died
    ok = died & np.isfinite(med_ho)
    if ok.sum() >= 10:
        err = np.abs(med_ho[ok] - th[ok])
        res["per_patient_point_estimate"] = {
            "n_evaluable_deaths": int(ok.sum()),
            "median_absolute_error_years": round(float(np.median(err)), 2),
            "mae_years": round(float(err.mean()), 2),
            "within_6_months": round(float((err <= 0.5).mean()), 3),
            "within_1_year": round(float((err <= 1.0).mean()), 3),
            "note": ("evaluated only on patients who died inside follow-up; a censored patient has no "
                     "known survival time to compare against, so including them would flatter the error"),
        }

    json.dump(res, open(OUT, "w"), indent=1)

    print("=== calibration on the sealed hold-out (%d patients, %d deaths) ==="
          % (res["n_holdout"], res["events_holdout"]))
    for h, c in cal.items():
        print("  %s horizon: overall predicted %.2f vs observed %.2f | mean |gap| across quartiles %s"
              % (h, c["overall_predicted"], c["overall_observed_km"], c["mean_abs_calibration_gap"]))
        for b in c["bins"]:
            print("      n=%-3d predicted %.2f  observed %.2f" % (b["n"], b["predicted"], b["observed_km"]))
    print("\n=== median survival by predicted risk tertile ===")
    for k, v in tert.items():
        print("  %-13s n=%-3d predicted %s y | observed %s y | error %s y"
              % (k, v["n"], v["predicted_median_years"], v["observed_km_median_years"], v["error_years"]))
    print("\n=== how long ONE patient lives, within a risk group (deaths only) ===")
    for k, v in spread.items():
        print("  %-13s median %.1f y, but 10th-90th percentile spans %.1f-%.1f y (n=%d)"
              % (k, v["median"], v["survival_years_p10_p90"][0], v["survival_years_p10_p90"][1], v["n_deaths"]))
    pp = res.get("per_patient_point_estimate")
    if pp:
        print("\n=== per-patient point estimate, among the %d who died ===" % pp["n_evaluable_deaths"])
        print("  median absolute error %.2f y | MAE %.2f y | within 6 months %.0f%% | within 1 year %.0f%%"
              % (pp["median_absolute_error_years"], pp["mae_years"],
                 100 * pp["within_6_months"], 100 * pp["within_1_year"]))
    print("\nwrote %s" % OUT)


if __name__ == "__main__":
    main()
