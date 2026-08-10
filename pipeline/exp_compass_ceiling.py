#!/usr/bin/env python3
"""How good could COMPASS possibly get? Estimate the ex-vivo assay's own reliability ceiling.

A model cannot predict a measurement more accurately than that measurement agrees with itself. Before
spending effort on architectures, measure the ceiling three ways:

  1. **Same-target drug pairs.** Two inhibitors of the same target, screened on the SAME specimens,
     measure (approximately) the same underlying biology with independent assay noise. Their observed
     correlation IS the reliability. A model that recovered the true biology perfectly would correlate
     with either measurement at only sqrt(reliability) — that is the ceiling.
  2. **Serial specimens from one patient.** How reproducible is a patient's drug-response profile at
     all? (Confounded by real biological change between timepoints, so this is a lower bound.)
  3. **Split-half over drugs within a family.** Correlate a specimen's mean response to a random half
     of a drug family against the other half — a direct internal-consistency estimate.

  python exp_compass_ceiling.py -> deliverables/exp_compass_ceiling.json
"""
import os, sys, json, time, warnings, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr

from amlmm.drug import data as D, targets as TG

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_compass_ceiling.json")
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")


def zmat(long):
    """specimen x drug matrix of within-drug robust z (higher = more resistant), double-centred so the
    patient main effect is removed — the same target the model is honestly judged against."""
    M = long.pivot_table(index="specimen", columns="inhibitor", values="auc", aggfunc="mean")
    Z = (M - M.median()) / (1.4826 * (M - M.median()).abs().median())
    return Z.sub(Z.median(axis=1), axis=0)          # remove the patient main effect


def main():
    t0 = time.time()
    d = np.load(BUNDLE, allow_pickle=True)
    ba = [str(s) for s in d["ba_samples"]]
    long = D.load(specimens=ba)
    elig = D.eligible_drugs(D.drug_summary(long))
    drugs = sorted(elig.loc[elig["tier"].isin(["primary", "wave_conditional"]), "inhibitor"])
    long = long[long["inhibitor"].isin(drugs)]
    Z = zmat(long)
    ann = TG.annotation()
    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "n_specimens": int(Z.shape[0]),
           "n_drugs": int(Z.shape[1])}

    # ---- 1. same-target pairs -------------------------------------------------------------
    pairs = []
    for a, b in itertools.combinations([c for c in Z.columns if c in ann], 2):
        ta, tb = set(ann[a]["targets"]), set(ann[b]["targets"])
        if not ta or not tb:
            continue
        jac = len(ta & tb) / len(ta | tb)
        if jac < 0.5:                                   # near-identical target profile only
            continue
        s = Z[[a, b]].dropna()
        if len(s) < 100:
            continue
        r = spearmanr(s[a], s[b]).statistic
        pairs.append({"a": a, "b": b, "jaccard": round(jac, 2), "n": int(len(s)),
                      "spearman": round(float(r), 4),
                      "family": ann[a]["family"]})
    pairs.sort(key=lambda p: -p["spearman"])
    rel = float(np.median([p["spearman"] for p in pairs])) if pairs else float("nan")
    res["same_target_pairs"] = {
        "n_pairs": len(pairs), "median_reliability_spearman": round(rel, 4),
        "implied_ceiling_spearman": round(float(np.sqrt(max(rel, 0))), 4),
        "top": pairs[:10], "bottom": pairs[-6:]}
    print("1. same-target drug pairs: %d pairs, median r = %.3f -> model ceiling ~%.3f Spearman"
          % (len(pairs), rel, np.sqrt(max(rel, 0))))
    for p in pairs[:5]:
        print("     %-26s vs %-26s r=%.3f (n=%d)" % (p["a"][:26], p["b"][:26], p["spearman"], p["n"]))

    # ---- 2. serial specimens from one patient ---------------------------------------------
    spec2sub = long.drop_duplicates("specimen").set_index("specimen")["subject"].to_dict()
    bysub = {}
    for s in Z.index:
        bysub.setdefault(spec2sub.get(s), []).append(s)
    rr = []
    for sub, ss in bysub.items():
        if len(ss) < 2:
            continue
        for a, b in itertools.combinations(ss, 2):
            v = Z.loc[[a, b]].T.dropna()
            if len(v) >= 60:
                rr.append(float(spearmanr(v.iloc[:, 0], v.iloc[:, 1]).statistic))
    res["serial_specimens"] = {"n_pairs": len(rr),
                               "median_profile_spearman": None if not rr else round(float(np.median(rr)), 4)}
    print("2. serial specimens from one patient: %d pairs, median profile r = %s"
          % (len(rr), res["serial_specimens"]["median_profile_spearman"]))

    # ---- 3. split-half within a drug family -----------------------------------------------
    fam = {}
    for dr in Z.columns:
        if dr in ann:
            fam.setdefault(ann[dr]["family_group"], []).append(dr)
    rng = np.random.RandomState(0)
    sh = []
    for f, ds in fam.items():
        if len(ds) < 6:
            continue
        vals = []
        for _ in range(30):
            p = rng.permutation(ds); h1, h2 = p[:len(p)//2], p[len(p)//2:]
            m1, m2 = Z[list(h1)].mean(1), Z[list(h2)].mean(1)
            v = pd.concat([m1, m2], axis=1).dropna()
            if len(v) >= 100:
                vals.append(float(spearmanr(v.iloc[:, 0], v.iloc[:, 1]).statistic))
        if vals:
            r_half = float(np.median(vals))
            # Spearman-Brown: reliability of the FULL family from a split half
            full = 2 * r_half / (1 + r_half) if r_half > -1 else np.nan
            sh.append({"family_group": f, "n_drugs": len(ds), "split_half_r": round(r_half, 4),
                       "spearman_brown_full": round(full, 4)})
    sh.sort(key=lambda x: -x["spearman_brown_full"])
    res["split_half_within_family"] = sh
    print("3. split-half within drug family (Spearman-Brown corrected):")
    for x in sh[:6]:
        print("     %-18s %2d drugs  half r=%.3f  full=%.3f" % (x["family_group"], x["n_drugs"],
                                                                 x["split_half_r"], x["spearman_brown_full"]))

    # ---- 4. how much of a drug's variance is shared with its family? ----------------------
    shared = []
    for dr in Z.columns:
        f = ann.get(dr, {}).get("family_group")
        sib = [x for x in fam.get(f, []) if x != dr]
        if len(sib) < 2:
            continue
        v = pd.concat([Z[dr], Z[sib].mean(1)], axis=1).dropna()
        if len(v) >= 100:
            shared.append({"inhibitor": dr, "family_group": f,
                           "r_with_family_mean": round(float(spearmanr(v.iloc[:, 0], v.iloc[:, 1]).statistic), 4)})
    shared.sort(key=lambda x: -x["r_with_family_mean"])
    res["drug_vs_family_mean"] = {"median": round(float(np.median([x["r_with_family_mean"] for x in shared])), 4),
                                  "top": shared[:8], "bottom": shared[-8:]}
    print("4. correlation of each drug with its family mean: median %.3f"
          % res["drug_vs_family_mean"]["median"])
    print("     most idiosyncratic (lowest shared signal):")
    for x in shared[-4:]:
        print("       %-28s %s  r=%.3f" % (x["inhibitor"][:28], x["family_group"], x["r_with_family_mean"]))

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
