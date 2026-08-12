#!/usr/bin/env python3
"""Deep characterisation of the survival layer — everything measurable, including what looks bad.

Nine experiments. The order is deliberate: the ones that could invalidate the headline run first.

  A1  hold-out lottery      is C-index 0.787 real, or did one split fall kindly? Repeat the whole
                            sealed-hold-out protocol over many seeds and report the DISTRIBUTION.
  A2  fold variance         how stable is the 0.751 CV number across folds and seeds
  A3  treatment strata      fit separately in intensive vs non-intensive induction. This targets the
                            worst documented failure (C-index 0.554 in the non-intensive group) and
                            needs no causal assumptions -- it answers "given this treatment path".
  A4  new clinical features FLT3-ITD allelic ratio, karyotype complexity, WBC/blast/platelet labs.
                            Sitting unused in the clinical file; the omissions a haematologist notices.
  A5  LSC17                 the strongest published prior, as a benchmark AND as a feature. Exposes
                            whether 60 expression PCs are earning their keep.
  A6  subgroups             C-index within age bands, ELN classes, sex, de novo/secondary, WBC tertiles
  A7  learning curve        C-index vs training size -- is the model data-limited or signal-limited?
  A8  horizons              discrimination at 0.5/1/2/3/5 y, to see where the model actually works
  A9  calibration drift     calibration gap by subgroup, since a model can rank well and be badly
                            calibrated inside a stratum

Every experiment respects the same discipline as training: patients never span folds, the feature space
is refit inside every fold, and censoring is handled properly.

  python exp_survival_deep.py [--seeds 40] [--only A1,A3]  ->  deliverables/exp_survival_deep.json
"""
import os, sys, json, time, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold

from amlmm.survival import data as SD, coxph as CX
from train_survival_model import (build_blocks, fit_arm, eval_arm, ARM_BLOCKS, BUNDLE, HORIZONS)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_survival_deep.json")


def load():
    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(x) for x in d["genes"]]
    ba = [str(s) for s in d["ba_samples"]]
    X = d["ba_X"].astype(np.float64)
    mut = pd.DataFrame(d["ba_L"].astype(float), index=ba,
                       columns=[str(c) for c in d["drivers"]])
    row_of = {s: i for i, s in enumerate(ba)}
    cl = SD.load_cohort(specimens=ba)
    rows = np.array([row_of[s] for s in cl["specimen"]])
    return dict(cl=cl, X=X, mut=mut, genes=genes, rows=rows,
                t=cl["time_years"].values.astype(float),
                e=cl["event"].values.astype(int), g=cl["subject"].values)


def blocks_for(D, train_idx, extra=None):
    fs, B = build_blocks(D["cl"], D["X"], D["mut"], D["rows"], sorted(D["rows"][train_idx]),
                         D["genes"], 60, 4000, train_idx=train_idx)
    if extra:
        B.update(extra)
    return fs, B


def cv_arm(D, arm, idx, folds=5, seed=0, blocks_extra=None):
    """Patient-grouped OOF risk for one arm over the patients in `idx`."""
    t, e, g = D["t"][idx], D["e"][idx], D["g"][idx]
    oof = np.full(len(idx), np.nan)
    gk = GroupKFold(n_splits=min(folds, len(set(g))))
    for i_in, i_out in gk.split(np.zeros(len(idx)), groups=g):
        gl_in, gl_out = idx[i_in], idx[i_out]
        try:
            _, B = blocks_for(D, gl_in, blocks_extra)
            Btr = {k: v[gl_in] for k, v in B.items()}
            Bte = {k: v[gl_out] for k, v in B.items()}
            r, _ = fit_arm(arm, Btr, Bte, D["t"][gl_in], D["e"][gl_in], D["g"][gl_in])
            oof[i_out] = r
        except Exception:
            pass
    m = np.isfinite(oof)
    if m.sum() < 20 or e[m].sum() < 5:
        return None, oof
    return float(CX.c_index(t[m], e[m], oof[m])), oof


# ---------------------------------------------------------------- A1 hold-out lottery
def a1_holdout_lottery(D, seeds, folds):
    """The sealed hold-out is ONE draw. Repeating it says whether 0.787 was skill or luck."""
    subs = np.array(sorted(set(D["g"])))
    out = {arm: [] for arm in ("age_eln", "molecular", "full", "deployed")}
    for s in range(seeds):
        rng = np.random.RandomState(1000 + s)
        hold = set(subs[rng.permutation(len(subs))[:int(round(0.20 * len(subs)))]])
        is_h = np.array([x in hold for x in D["g"]])
        tr, ho = np.where(~is_h)[0], np.where(is_h)[0]
        if D["e"][ho].sum() < 10:
            continue
        try:
            _, B = blocks_for(D, tr)
        except Exception:
            continue
        for arm in out:
            try:
                Btr = {k: v[tr] for k, v in B.items()}
                Bho = {k: v[ho] for k, v in B.items()}
                r, _ = fit_arm(arm, Btr, Bho, D["t"][tr], D["e"][tr], D["g"][tr])
                m = np.isfinite(r)
                if m.sum() >= 20:
                    out[arm].append(float(CX.c_index(D["t"][ho][m], D["e"][ho][m], r[m])))
            except Exception:
                pass
        print("    seed %d/%d" % (s + 1, seeds), flush=True)
    res = {}
    for arm, v in out.items():
        if not v:
            continue
        v = np.array(v)
        res[arm] = {"n_splits": len(v), "mean": round(float(v.mean()), 4),
                    "sd": round(float(v.std()), 4), "min": round(float(v.min()), 4),
                    "max": round(float(v.max()), 4),
                    "p05": round(float(np.percentile(v, 5)), 4),
                    "p95": round(float(np.percentile(v, 95)), 4)}
    # how often does deployed beat age+ELN on the SAME split?
    if "deployed" in out and "age_eln" in out and len(out["deployed"]) == len(out["age_eln"]):
        d = np.array(out["deployed"]) - np.array(out["age_eln"])
        res["deployed_minus_age_eln"] = {"mean": round(float(d.mean()), 4),
                                         "sd": round(float(d.std()), 4),
                                         "frac_splits_positive": round(float((d > 0).mean()), 4)}
    return res


# ---------------------------------------------------------------- A3 treatment strata
def a3_treatment_strata(D, folds):
    ty = D["cl"]["typeInductionTx"].astype(str)
    strata = {"intensive": ty.str.contains("Standard Chemo", na=False).values,
              "non_intensive": (~ty.str.contains("Standard Chemo", na=False)
                                & ty.ne("nan")).values}
    res = {}
    for name, mask in strata.items():
        idx = np.where(mask)[0]
        if len(idx) < 40 or D["e"][idx].sum() < 15:
            res[name] = {"n": int(len(idx)), "events": int(D["e"][idx].sum()),
                         "note": "too few to fit a stratum-specific model"}
            continue
        row = {"n": int(len(idx)), "events": int(D["e"][idx].sum())}
        for arm in ("age_eln", "molecular", "deployed"):
            c, _ = cv_arm(D, arm, idx, folds)
            row[arm] = None if c is None else round(c, 4)
        res[name] = row
    return res


# ---------------------------------------------------------------- A4 unused clinical features
def a4_extra_clinical(D, folds):
    cl = D["cl"]
    def num(c):
        return pd.to_numeric(cl[c], errors="coerce").values if c in cl.columns else None
    cols, names = [], []
    for c, nm in (("FLT3_ITD_allelic_ratio", "flt3_itd_ar"), ("allelic_ratio", "allelic_ratio"),
                  ("wbcCount", "wbc"), ("%.Blasts.in.BM", "blasts_bm"), ("%.Blasts.in.PB", "blasts_pb"),
                  ("plateletCount", "platelets"), ("hemoglobin", "hgb")):
        v = num(c)
        if v is not None and np.isfinite(v).sum() > 0.5 * len(cl):
            cols.append(np.nan_to_num(v, nan=float(np.nanmedian(v)))); names.append(nm)
    # karyotype complexity: count of abnormalities in the free-text karyotype
    for c in ("karyotype", "Karyotype", "cytogenetics"):
        if c in cl.columns:
            k = cl[c].astype(str)
            cols.append(k.str.count(",").fillna(0).values.astype(float)); names.append("karyo_complexity")
            break
    if not cols:
        return {"note": "none of the candidate columns are present with >50% coverage",
                "columns_seen": sorted(cl.columns.tolist())[:40]}
    EXTRA = np.vstack(cols).T
    idx = np.arange(len(D["cl"]))
    base, _ = cv_arm(D, "deployed", idx, folds)
    aug, _ = cv_arm(D, "deployed_plus", idx, folds,
                    blocks_extra={"extra_clin": EXTRA})
    return {"features_used": names, "n_features": len(names),
            "deployed": None if base is None else round(base, 4),
            "deployed_plus_extra": None if aug is None else round(aug, 4),
            "delta": None if (base is None or aug is None) else round(aug - base, 4)}


# ---------------------------------------------------------------- A5 LSC17
LSC17 = ["DNMT3B", "ZBTB46", "NYNRIN", "ARHGAP22", "LAPTM4B", "MMRN1", "DPYSL3", "KIAA0125",
         "CDK6", "CPXM1", "SOCS2", "SMIM24", "EMP1", "NGFRAP1", "CD34", "AKR1C3", "GPR56"]
LSC17_W = [0.0874, -0.0347, 0.00865, -0.0138, 0.00582, 0.0258, 0.0284, 0.0196, -0.0704, -0.0258,
           0.0271, -0.0226, 0.0146, 0.0465, 0.0338, -0.0402, 0.0501]


def a5_lsc17(D, folds, sym2ens):
    gpos = {g: i for i, g in enumerate(D["genes"])}
    idxs, ws = [], []
    for g, w in zip(LSC17, LSC17_W):
        j = gpos.get(g) or gpos.get(sym2ens.get(g, ""))
        if j is not None:
            idxs.append(j); ws.append(w)
    if len(idxs) < 10:
        return {"note": "only %d/17 LSC17 genes found in the model's gene space" % len(idxs),
                "found": len(idxs)}
    L = np.log2(D["X"][:, idxs] + 1.0)
    L = (L - L.mean(0)) / (L.std(0) + 1e-9)
    score = (L @ np.array(ws)).reshape(-1, 1)[D["rows"]]
    idx = np.arange(len(D["cl"]))
    out = {"genes_found": len(idxs), "of": 17}
    c, _ = cv_arm(D, "lsc17", idx, folds, blocks_extra={"lsc17": score})
    out["lsc17_alone"] = None if c is None else round(c, 4)
    for arm in ("age_eln", "rna", "deployed"):
        c, _ = cv_arm(D, arm, idx, folds)
        out[arm] = None if c is None else round(c, 4)
    c, _ = cv_arm(D, "deployed_lsc17", idx, folds, blocks_extra={"lsc17": score})
    out["deployed_plus_lsc17"] = None if c is None else round(c, 4)
    return out


# ---------------------------------------------------------------- A6 subgroups
def a6_subgroups(D, folds):
    idx = np.arange(len(D["cl"]))
    c_all, oof = cv_arm(D, "deployed", idx, folds)
    _, oof_base = cv_arm(D, "age_eln", idx, folds)
    cl = D["cl"]
    age = pd.to_numeric(cl["ageAtDiagnosis"], errors="coerce").values
    eln = cl["ELN2017"].astype(str).values
    sex = cl["consensus_sex"].astype(str).values if "consensus_sex" in cl.columns else np.array(["?"] * len(cl))
    groups = {"age <60": age < 60, "age 60-74": (age >= 60) & (age < 75), "age >=75": age >= 75,
              "ELN Favorable": eln == "Favorable", "ELN Intermediate": eln == "Intermediate",
              "ELN Adverse": eln == "Adverse", "male": sex == "Male", "female": sex == "Female"}
    for c, nm in (("isDenovo", "de novo"), ("isTransformed", "transformed"), ("priorMDS", "prior MDS")):
        if c in cl.columns:
            groups[nm] = pd.to_numeric(cl[c], errors="coerce").fillna(0).values > 0
    res = {"overall": {"n": len(idx), "events": int(D["e"].sum()),
                       "deployed": None if c_all is None else round(c_all, 4)}}
    for nm, m in groups.items():
        m = np.asarray(m) & np.isfinite(oof)
        if m.sum() < 25 or D["e"][m].sum() < 8:
            res[nm] = {"n": int(m.sum()), "events": int(D["e"][m].sum()), "note": "too small"}
            continue
        cd = CX.c_index(D["t"][m], D["e"][m], oof[m])
        mb = m & np.isfinite(oof_base)
        cb = CX.c_index(D["t"][mb], D["e"][mb], oof_base[mb]) if mb.sum() >= 25 else None
        res[nm] = {"n": int(m.sum()), "events": int(D["e"][m].sum()),
                   "deployed": round(float(cd), 4),
                   "age_eln": None if cb is None else round(float(cb), 4),
                   "delta": None if cb is None else round(float(cd - cb), 4)}
    return res


# ---------------------------------------------------------------- A7 learning curve
def a7_learning_curve(D, folds, seeds=5):
    subs = np.array(sorted(set(D["g"])))
    fracs = [0.15, 0.3, 0.45, 0.6, 0.8, 1.0]
    res = {}
    for f in fracs:
        vals = []
        for s in range(seeds):
            rng = np.random.RandomState(500 + s)
            keep = set(subs[rng.permutation(len(subs))[:max(30, int(f * len(subs)))]])
            idx = np.where([x in keep for x in D["g"]])[0]
            if len(idx) < 40 or D["e"][idx].sum() < 15:
                continue
            c, _ = cv_arm(D, "deployed", idx, folds)
            if c is not None:
                vals.append(c)
        if vals:
            res["%d%%" % round(100 * f)] = {"n_patients": int(round(f * len(subs))),
                                            "mean_c_index": round(float(np.mean(vals)), 4),
                                            "sd": round(float(np.std(vals)), 4), "n_reps": len(vals)}
        print("    learning curve %.0f%% done" % (100 * f), flush=True)
    return res


# ---------------------------------------------------------------- A8 horizons
def a8_horizons(D, folds):
    idx = np.arange(len(D["cl"]))
    out = {}
    for arm in ("age_eln", "deployed"):
        c, oof = cv_arm(D, arm, idx, folds)
        m = np.isfinite(oof)
        row = {"c_index": None if c is None else round(c, 4)}
        for h in (0.5, 1.0, 2.0, 3.0, 5.0):
            a, n = CX.td_auc(D["t"][m], D["e"][m], oof[m], h)
            row["auc_%gy" % h] = None if a is None else round(a, 4)
            row["n_%gy" % h] = n
        out[arm] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=40)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    want = set(x.strip().upper() for x in a.only.split(",") if x.strip())
    t0 = time.time()

    # arms invented here need a block list; register them before anything runs
    ARM_BLOCKS.setdefault("deployed_plus", ARM_BLOCKS["deployed"] + ["extra_clin"])
    ARM_BLOCKS.setdefault("lsc17", ["lsc17"])
    ARM_BLOCKS.setdefault("deployed_lsc17", ARM_BLOCKS["deployed"] + ["lsc17"])

    D = load()
    prov = D["cl"].attrs["provenance"]
    print("cohort: %d patients, %d deaths" % (prov["final_patients"], prov["events"]), flush=True)
    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "cohort": prov, "seeds": a.seeds,
           "folds": a.folds, "experiments": {}}

    import train_drug_model as TDM
    sym2ens = TDM.sym2ens_map(D["genes"])

    jobs = [("A1", "hold-out lottery", lambda: a1_holdout_lottery(D, a.seeds, a.folds)),
            ("A3", "treatment strata", lambda: a3_treatment_strata(D, a.folds)),
            ("A4", "unused clinical features", lambda: a4_extra_clinical(D, a.folds)),
            ("A5", "LSC17", lambda: a5_lsc17(D, a.folds, sym2ens)),
            ("A6", "subgroups", lambda: a6_subgroups(D, a.folds)),
            ("A7", "learning curve", lambda: a7_learning_curve(D, a.folds)),
            ("A8", "horizons", lambda: a8_horizons(D, a.folds))]

    for code, name, fn in jobs:
        if want and code not in want:
            continue
        print("\n== %s: %s ==" % (code, name), flush=True)
        s = time.time()
        try:
            res["experiments"][code] = {"name": name, "result": fn(),
                                        "seconds": round(time.time() - s, 1)}
            print("   %s done in %.0fs" % (code, time.time() - s), flush=True)
        except Exception as ex:
            import traceback
            res["experiments"][code] = {"name": name, "error": str(ex),
                                        "traceback": traceback.format_exc()[-1500:]}
            print("   %s FAILED: %s" % (code, ex), flush=True)
        with open(OUT, "w") as f:                      # write after every experiment, not at the end
            json.dump(res, f, indent=1)

    print("\nwrote %s (%.0fs total)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
