#!/usr/bin/env python3
"""Validate Model B (state-resolved response) on the single-cell atlas.

There is no ground truth for state-level ex-vivo response -- BeatAML measured bulk mononuclear cells --
so this cannot be a supervised accuracy check, and pretending otherwise would be the dishonest move.
What CAN be tested, and is tested here:

  1. **A pre-registered biological expectation.** Venetoclax response in AML is repeatedly reported to
     track differentiation state: monocytic AML is relatively resistant (MCL1-dependent), primitive /
     stem-like AML relatively sensitive (BCL2-dependent). Model B never saw a cell-state label during
     training -- it was fitted entirely on BeatAML bulk. If it reproduces that contrast when handed
     cell states it has never been told about, that is a real, falsifiable result. If it does not, we
     say so.
  2. **A negative control.** The same contrast is computed for every other modelled drug. If primitive
     states simply score as "sensitive to everything", the venetoclax result means nothing -- so the
     venetoclax effect is reported as a rank among all inhibitors, not in isolation.
  3. **Internal consistency.** sens_bulk vs the abundance-weighted mean of the state predictions.
  4. **Out-of-distribution distance.** How far the single-cell bulk-equivalents sit from the BeatAML
     training distribution the model was fitted on.

  python validate_state_response.py [--n-samples 60]  ->  deliverables/state_response_validation.json
"""
import os, sys, json, time, argparse, warnings, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
from scipy.stats import wilcoxon, spearmanr

from amlmm.drug import statemodel as SM, targets as TG

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
H5 = os.path.join(ROOT, "data", "RNA", "pseudobulk_counts_hashed.h5ad")
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
MODEL = os.path.join(HERE, "drug_response_model.pkl")
OUT = os.path.join(ROOT, "deliverables", "state_response_validation.json")
STATE_COL = "Hs-BM-titrated-reference-centroid"


def read_atlas(samples=None):
    """Row-subset read: the obs columns are small, so pick the wanted rows first and pull only those
    out of the CSR. Materialising the whole 1 GB matrix to score 80 samples is the difference between
    a minute and half an hour."""
    from amlmm.drug.h5rows import obs_column, var_index, read_rows
    with h5py.File(H5, "r") as f:
        ds, sm, st = obs_column(f, "Dataset"), obs_column(f, "Sample"), obs_column(f, STATE_COL)
        nc = np.asarray(obs_column(f, "n_cells"), dtype=float)
        genes = var_index(f)
    key = np.array(["%s::%s" % (a, b) for a, b in zip(ds, sm)])
    rows = np.arange(len(key)) if samples is None else np.where(np.isin(key, list(samples)))[0]
    _, X = read_rows(H5, rows)
    return X, genes, key[rows], st[rows], nc[rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-samples", type=int, default=60)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()

    with open(MODEL, "rb") as f:
        mod = pickle.load(f)
    sr = SM.StateResponse(mod)
    drugs = list(mod.drug_models)

    import h5py as _h5
    from amlmm.drug.h5rows import obs_column as _oc
    with _h5.File(H5, "r") as _f:
        _all = sorted(set("%s::%s" % (a_, b_) for a_, b_ in zip(_oc(_f, "Dataset"), _oc(_f, "Sample"))))
    rng = np.random.RandomState(a.seed)
    pick = list(rng.permutation(_all))[:a.n_samples]
    X, genes, key, st, nc = read_atlas(set(pick))
    print("atlas: %d samples total, evaluating %d (%s rows read) | %d drugs"
          % (len(_all), len(pick), X.shape, len(drugs)))

    prim_minus_mono = {d: [] for d in drugs}
    consistency, ood, per_sample = [], [], {}
    for si, s in enumerate(pick):
        m = key == s
        counts = pd.DataFrame(X[m], index=st[m], columns=genes).groupby(level=0).sum()
        ncells = pd.Series(nc[m], index=st[m]).groupby(level=0).sum()
        r = sr.predict(counts, ncells, drugs=drugs)
        if not r["per_drug"]:
            continue
        info = pd.DataFrame(r["states"]).set_index("state")
        prim = info.index[info["group"].isin(SM.LSC_GROUPS)]
        mono = info.index[info["group"] == "monocytic"]
        if len(prim) >= 2 and len(mono) >= 2:
            for d in drugs:
                ps = r["per_state"].get(d)
                if not ps:
                    continue
                p = np.mean([ps[x] for x in prim if x in ps])
                mo = np.mean([ps[x] for x in mono if x in ps])
                if np.isfinite(p) and np.isfinite(mo):
                    prim_minus_mono[d].append(float(p - mo))
        pdd = r["per_drug"]
        consistency += [(v["sens_bulk"], v["sens_weighted"]) for v in pdd.values()]
        per_sample[s] = {"n_states": len(info), "n_cells": int(info["n_cells"].sum()),
                         "venetoclax": pdd.get("Venetoclax")}
        if si and si % 20 == 0:
            print("   %d/%d (%.0fs)" % (si, len(pick), time.time() - t0))

    # ---- 1 & 2: the primitive-minus-monocytic contrast, venetoclax vs every other inhibitor ----
    rows = []
    for d, v in prim_minus_mono.items():
        if len(v) < 10:
            continue
        v = np.array(v)
        try:
            p = float(wilcoxon(v).pvalue)
        except Exception:
            p = None
        rows.append({"inhibitor": d, "n_samples": int(len(v)),
                     "mean_primitive_minus_monocytic": round(float(v.mean()), 4),
                     "frac_positive": round(float((v > 0).mean()), 3), "wilcoxon_p": p,
                     "clinical_tier": TG.get(d)["clinical_tier"], "family_group": TG.get(d)["family_group"]})
    contrast = pd.DataFrame(rows).sort_values("mean_primitive_minus_monocytic", ascending=False)
    contrast["rank"] = np.arange(1, len(contrast) + 1)
    ven = contrast[contrast["inhibitor"] == "Venetoclax"]

    cons = np.array(consistency)
    r_cons = spearmanr(cons[:, 0], cons[:, 1]) if len(cons) > 10 else None

    # ---- 4: OOD distance of the atlas bulk-equivalents from BeatAML training space ----
    d = np.load(BUNDLE, allow_pickle=True)
    fs = mod.fs
    Zba = fs.z(d["ba_X"].astype(float), "beataml")
    Zsc = fs.z(d["sc_X"].astype(float), "beataml")
    Pba, Psc = fs.pca.transform(Zba[:, fs.sel]), fs.pca.transform(Zsc[:, fs.sel])
    mu, sd = Pba.mean(0), Pba.std(0); sd[sd == 0] = 1.0
    dist_ba = np.sqrt((((Pba - mu) / sd) ** 2).mean(1))
    dist_sc = np.sqrt((((Psc - mu) / sd) ** 2).mean(1))

    res = {
        "generated": time.strftime("%Y-%m-%d %H:%M"),
        "n_samples_evaluated": len(per_sample), "n_drugs": len(drugs),
        "primitive_vs_monocytic": {
            "definition": "mean predicted sensitivity in HSC/MPP+LMPP/CLP+GMP states minus monocytic states",
            "venetoclax": (ven.iloc[0].to_dict() if len(ven) else None),
            "venetoclax_rank_of": int(len(contrast)),
            "top10": contrast.head(10).to_dict("records"),
            "bottom10": contrast.tail(10).to_dict("records"),
        },
        "internal_consistency": {
            "n_pairs": int(len(cons)),
            "spearman_bulk_vs_weighted": None if r_cons is None else round(float(r_cons.statistic), 4),
            "mean_abs_difference": round(float(np.mean(np.abs(cons[:, 0] - cons[:, 1]))), 4),
        },
        "out_of_distribution": {
            "metric": "RMS z-distance in BeatAML PCA space (per-PC standardised)",
            "beataml_median": round(float(np.median(dist_ba)), 3),
            "beataml_p95": round(float(np.percentile(dist_ba, 95)), 3),
            "singlecell_median": round(float(np.median(dist_sc)), 3),
            "singlecell_frac_beyond_beataml_p95":
                round(float((dist_sc > np.percentile(dist_ba, 95)).mean()), 3),
        },
    }
    json.dump(res, open(OUT, "w"), indent=1)

    print("\n== primitive minus monocytic predicted sensitivity ==")
    print(contrast.head(12)[["rank", "inhibitor", "mean_primitive_minus_monocytic",
                             "frac_positive", "wilcoxon_p", "clinical_tier"]].to_string(index=False))
    if len(ven):
        v = ven.iloc[0]
        print("\nVENETOCLAX: rank %d/%d, mean primitive-minus-monocytic %+.3f, %.0f%% of samples positive, p=%s"
              % (v["rank"], len(contrast), v["mean_primitive_minus_monocytic"],
                 100 * v["frac_positive"], v["wilcoxon_p"]))
    print("\ninternal consistency  Spearman(bulk, weighted) = %s | mean |diff| %.3f"
          % (res["internal_consistency"]["spearman_bulk_vs_weighted"],
             res["internal_consistency"]["mean_abs_difference"]))
    print("OOD: single-cell median distance %.2f vs BeatAML %.2f; %.0f%% beyond BeatAML p95"
          % (res["out_of_distribution"]["singlecell_median"], res["out_of_distribution"]["beataml_median"],
             100 * res["out_of_distribution"]["singlecell_frac_beyond_beataml_p95"]))
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
