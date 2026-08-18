#!/usr/bin/env python3
"""Train and honestly evaluate the MOSAIC-AML survival layer on BeatAML2.

The question is not "what C-index can we get". Age and ELN 2017 already predict AML survival, are free,
and are what a clinician has on day one. So the number that matters is **the gain over age + ELN**, and
every arm below exists to make that comparison unavoidable:

    age            age alone
    eln            ELN 2017 alone
    age_eln        THE BASELINE TO BEAT
    clin           the full clinical block (age, ELN, labs, prior disease)
    rna/state/mut  each molecular block on its own
    molecular      the three molecular blocks fused, no clinical information at all
    full           everything fused

Discipline is the same as the rest of the platform: patients never span folds, a sealed hold-out of
patients is drawn once and never touched, and the feature space (PCA, imputation medians, z-reference)
is refit inside every fold. Censoring is handled properly throughout -- a patient last known alive
contributes as censored at their follow-up time, never as a death.

  python train_survival_model.py  ->  pipeline/survival_model.pkl
                                      deliverables/survival_model_card.json
"""
import os, sys, json, time, pickle, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold

from amlmm.survival import data as SD, coxph as CX
from amlmm.survival.model import SurvivalModel
from amlmm.drug.features import FeatureSpace

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
OUTP = os.path.join(HERE, "survival_model.pkl")
CARD = os.path.join(ROOT, "deliverables", "survival_model_card.json")

HORIZONS = [1.0, 2.0, 5.0]                       # years
# Deployed after the sweep measured them: baseline induction type takes C-index 0.726 -> 0.750 (the
# single biggest gain anywhere in the platform) and an age spline adds a further ~0.007. Only
# DIAGNOSIS-TIME treatment information is used -- n_regimens correlates 0.571 with follow-up time
# because a patient must survive to accumulate regimens, and including it would read the outcome.
ARMS = ["age", "eln", "age_eln", "clin", "rna", "state", "mut", "molecular", "full", "deployed"]
ARM_BLOCKS = {"age": ["age"], "eln": ["eln"], "age_eln": ["age_eln"], "clin": ["clin"],
              "rna": ["rna"], "state": ["state"], "mut": ["mut"],
              "molecular": ["rna", "state", "mut"], "full": ["rna", "state", "mut", "clin"],
              "deployed": ["rna", "state", "mut", "clin", "age_spline", "txbase"]}


def build_blocks(cl, X_lin, mut_all, rows, train_rows, genes, n_pc, var_genes, train_idx=None):
    """All feature blocks for this cohort, with the expression space fit on `train_rows` only."""
    fs = FeatureSpace(genes, n_pc=n_pc, var_genes=var_genes)
    fs.add_reference("beataml", X_lin[train_rows])
    Z = fs.z(X_lin, "beataml")
    fs.fit(Z[train_rows], None, None)
    P = fs.pca.transform(Z[:, fs.sel])[rows]
    S, _ = fs._state_block(Z)
    S = S[rows]
    # ba_L is 14% NaN ("not assayed"), every row has at least one, and one column is entirely NaN.
    # Feeding that to Cox fails silently and — before this was caught — took every FUSED arm down with
    # it, because the fusion's "rows all blocks scored" mask became empty. Drop the dead and the too-rare
    # columns, zero-fill the rest (not-assayed reads as not-observed, as elsewhere in the platform), and
    # keep one column recording how much was missing so the model can discount a poorly-assayed sample.
    Mraw = mut_all.iloc[rows].values.astype(float)
    keep = (~np.isnan(Mraw).all(0)) & (np.nansum(Mraw, 0) >= 10)
    miss = np.isnan(Mraw[:, keep]).mean(1, keepdims=True)
    M = np.hstack([np.nan_to_num(Mraw[:, keep], nan=0.0), miss])
    build_blocks.mut_columns = list(np.asarray(mut_all.columns)[keep]) + ["__frac_missing"]
    C, _ = SD.clinical_block(cl)
    AE, _ = SD.age_eln_block(cl)
    age = AE[:, [0]]
    eln = AE[:, 1:]
    # age spline (piecewise-cubic at the training quartiles) and DIAGNOSIS-TIME treatment only
    kn = np.quantile(age[train_idx] if train_idx is not None else age, [0.25, 0.5, 0.75])
    # the knots are part of the fitted model, not a constant: coefficients fitted against knots at the
    # training quartiles are wrong if inference rebuilds the basis somewhere else. Recorded here and
    # carried in the bundle so survival_layer.py uses the real ones instead of a hardcoded guess.
    build_blocks.age_knots = [float(q) for q in kn]
    AGS = np.hstack([age] + [np.clip(age - q, 0, None) ** 3 for q in kn])
    ty = cl["typeInductionTx"].astype(str)
    TXB = np.vstack([ty.str.contains("Standard Chemo").astype(float).values,
                     ty.str.contains("Palliative").astype(float).values,
                     ty.eq("nan").astype(float).values]).T
    return fs, {"rna": P, "state": S, "mut": M, "clin": C, "age": age, "eln": eln, "age_eln": AE,
                "age_spline": AGS, "txbase": TXB}


def fit_arm(arm, blocks_tr, blocks_te, t, e, g):
    """One arm: single-block Cox, or the stacked SurvivalModel when the arm has several blocks."""
    names = ARM_BLOCKS[arm]
    if len(names) == 1:
        b = names[0]
        m = CX.CoxPH(alpha=0.30 if b == "rna" else 0.10).fit(blocks_tr[b], t, e)
        return m.risk(blocks_te[b]), m
    sm = SurvivalModel().fit({b: blocks_tr[b] for b in names}, t, e, g)
    return sm.risk({b: blocks_te[b] for b in names}), sm


def eval_arm(t, e, risk, horizons=HORIZONS, surv=None):
    # An arm that failed in one fold leaves NaNs behind; score it on the patients it did cover and say
    # so, rather than crashing or silently imputing.
    m = np.isfinite(risk)
    if m.sum() < 20:
        return {"c_index": None, "n_scored": int(m.sum()), "note": "too few scored patients"}
    t, e, risk = np.asarray(t)[m], np.asarray(e)[m], np.asarray(risk)[m]
    if surv is not None:
        surv = np.asarray(surv)[m]
    out = {"c_index": round(CX.c_index(t, e, risk), 4), "n": int(len(t)), "events": int(e.sum()),
           "n_scored": int(m.sum()), "coverage": round(float(m.mean()), 3)}
    for h in horizons:
        a, n = CX.td_auc(t, e, risk, h)
        out["auc_%gy" % h] = None if a is None else round(a, 4)
        out["auc_%gy_n" % h] = n
        if surv is not None:
            out["brier_%gy" % h] = round(CX.ipcw_brier(t, e, surv[:, horizons.index(h)], h), 4)
    return out


def boot_delta_c(t, e, r_new, r_base, B=1000, seed=0):
    """Bootstrap the C-index DIFFERENCE over patients — the honest way to ask whether the molecular
    model adds anything, since the two arms are scored on the very same patients."""
    rng = np.random.RandomState(seed)
    m = np.isfinite(r_new) & np.isfinite(r_base)          # compare on the patients BOTH arms scored
    t, e, r_new, r_base = np.asarray(t)[m], np.asarray(e)[m], np.asarray(r_new)[m], np.asarray(r_base)[m]
    n = len(t)
    d = []
    for _ in range(B):
        i = rng.randint(0, n, n)
        if e[i].sum() < 5:
            continue
        d.append(CX.c_index(t[i], e[i], r_new[i]) - CX.c_index(t[i], e[i], r_base[i]))
    d = np.asarray(d)
    return {"delta_c": round(float(CX.c_index(t, e, r_new) - CX.c_index(t, e, r_base)), 4),
            "ci95": [round(float(np.percentile(d, 2.5)), 4), round(float(np.percentile(d, 97.5)), 4)],
            "p_gt0": round(float((d <= 0).mean()), 4), "B": int(len(d))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--holdout", type=float, default=0.20)
    ap.add_argument("--n-pc", type=int, default=60)
    ap.add_argument("--var-genes", type=int, default=4000)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()

    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(x) for x in d["genes"]]
    ba = [str(s) for s in d["ba_samples"]]
    X_lin = d["ba_X"].astype(np.float64)
    mut_all = pd.DataFrame(d["ba_L"].astype(float), index=ba, columns=[str(c) for c in d["drivers"]])
    row_of = {s: i for i, s in enumerate(ba)}

    cl = SD.load_cohort(specimens=ba)
    prov = cl.attrs["provenance"]
    rows = np.array([row_of[s] for s in cl["specimen"]])
    t = cl["time_years"].values.astype(float)
    e = cl["event"].values.astype(int)
    g = cl["subject"].values

    print("MOSAIC-AML survival layer")
    print("  cohort: %d patients, %d deaths (%.0f%%), %d censored; median follow-up of the living %.1f y"
          % (prov["final_patients"], prov["events"], 100 * prov["events"] / prov["final_patients"],
             prov["censored"], prov["median_followup_years_censored"]))
    print("  dropped: %d unknown vital status; non-initial-diagnosis specimens excluded"
          % prov["dropped_unknown_vital_status"])

    # ---------------- sealed hold-out of PATIENTS ----------------
    subs = np.array(sorted(set(g)))
    rng = np.random.RandomState(a.seed)
    perm = rng.permutation(len(subs))
    hold = set(subs[perm[:int(round(a.holdout * len(subs)))]])
    is_hold = np.array([s in hold for s in g])
    print("  sealed hold-out: %d patients (%.0f%%), %d deaths"
          % (is_hold.sum(), 100 * is_hold.mean(), int(e[is_hold].sum())))

    tr = ~is_hold
    idx_tr = np.where(tr)[0]

    # ---------------- patient-grouped CV inside the training partition ----------------
    oof = {arm: np.full(tr.sum(), np.nan) for arm in ARMS}
    gk = GroupKFold(n_splits=a.folds)
    for k, (i_in, i_out) in enumerate(gk.split(np.zeros(tr.sum()), groups=g[tr])):
        gl_in, gl_out = idx_tr[i_in], idx_tr[i_out]
        _, B = build_blocks(cl, X_lin, mut_all, rows, sorted(rows[gl_in]), genes, a.n_pc, a.var_genes,
                            train_idx=gl_in)
        Btr = {kk: v[gl_in] for kk, v in B.items()}
        Bte = {kk: v[gl_out] for kk, v in B.items()}
        for arm in ARMS:
            try:
                r, _ = fit_arm(arm, Btr, Bte, t[gl_in], e[gl_in], g[gl_in])
                oof[arm][i_out] = r
            except Exception as ex:
                print("   fold %d arm %s failed: %s" % (k + 1, arm, ex))
        print("   fold %d/%d done (%.0fs)" % (k + 1, a.folds, time.time() - t0))

    t_tr, e_tr = t[tr], e[tr]
    res = {"cohort": prov, "horizons": HORIZONS, "folds": a.folds,
           "holdout_patients": int(is_hold.sum()), "holdout_events": int(e[is_hold].sum()),
           "oof": {}, "holdout": {}, "incremental": {}}
    for arm in ARMS:
        res["oof"][arm] = eval_arm(t_tr, e_tr, oof[arm])

    print("\n== donor-grouped CV inside the training partition (n=%d, %d deaths) =="
          % (tr.sum(), int(e_tr.sum())))
    print("  %-11s %8s %8s %8s %8s" % ("arm", "C-index", "AUC 1y", "AUC 2y", "AUC 5y"))
    for arm in ARMS:
        o = res["oof"][arm]
        if o.get("c_index") is None:
            print("  %-11s   (not scored: %s)" % (arm, o.get("note"))); continue
        print("  %-11s %8.3f %8s %8s %8s   cov %.2f" % (arm, o["c_index"], o["auc_1y"], o["auc_2y"],
                                                        o["auc_5y"], o["coverage"]))

    # ---------------- incremental value over age + ELN ----------------
    base = oof["age_eln"]
    # `deployed` belongs here: it is the arm inference actually selects when a patient has clinical
    # covariates, and without it every such report showed a blank where its gain over age + ELN goes.
    for arm in ("clin", "rna", "state", "mut", "molecular", "full", "deployed"):
        res["incremental"][arm] = boot_delta_c(t_tr, e_tr, oof[arm], base, B=a.boot, seed=a.seed)
    print("\n== incremental C-index over age + ELN 2017 (bootstrap over patients) ==")
    for arm, v in res["incremental"].items():
        print("  %-11s %+.4f  95%% CI [%+.4f, %+.4f]   P(no gain) = %.3f"
              % (arm, v["delta_c"], v["ci95"][0], v["ci95"][1], v["p_gt0"]))

    # ---------------- deployed model + sealed hold-out ----------------
    fs, B = build_blocks(cl, X_lin, mut_all, rows, sorted(rows[tr]), genes, a.n_pc, a.var_genes,
                         train_idx=idx_tr)
    Btr = {k: v[tr] for k, v in B.items()}
    Bho = {k: v[is_hold] for k, v in B.items()}
    fitted = {}
    for arm in ARMS:
        try:
            r, m = fit_arm(arm, Btr, Bho, t[tr], e[tr], g[tr])
            fitted[arm] = m
            sm = m.survival({b: Bho[b] for b in ARM_BLOCKS[arm]}, HORIZONS) if isinstance(m, SurvivalModel) \
                else m.survival(Bho[ARM_BLOCKS[arm][0]], HORIZONS)
            res["holdout"][arm] = eval_arm(t[is_hold], e[is_hold], r, surv=sm)
        except Exception as ex:
            print("   holdout arm %s failed: %s" % (arm, ex))
    print("\n== sealed hold-out (%d patients, %d deaths) ==" % (is_hold.sum(), int(e[is_hold].sum())))
    print("  %-11s %8s %8s %8s" % ("arm", "C-index", "AUC 2y", "Brier 2y"))
    for arm in ARMS:
        h = res["holdout"].get(arm)
        if h and h.get("c_index") is not None:
            print("  %-11s %8.3f %8s %8s" % (arm, h["c_index"], h["auc_2y"], h.get("brier_2y")))

    # ---------------- risk-group separation ----------------
    best = "deployed" if "deployed" in fitted else ("full" if "full" in fitted else "molecular")
    r_ho = fitted[best].risk({b: Bho[b] for b in ARM_BLOCKS[best]}) if isinstance(fitted[best], SurvivalModel) \
        else fitted[best].risk(Bho[ARM_BLOCKS[best][0]])
    cut = np.quantile(oof[best], [1 / 3, 2 / 3])
    grp = np.digitize(r_ho, cut)
    km_curves = {}
    for gi, name in enumerate(("low", "intermediate", "high")):
        m = grp == gi
        if m.sum() >= 5:
            ts, ss = CX.km(t[is_hold][m], e[is_hold][m])
            km_curves[name] = {"n": int(m.sum()), "events": int(e[is_hold][m].sum()),
                               "times": [round(float(x), 3) for x in ts],
                               "survival": [round(float(x), 4) for x in ss]}
    lo_hi = (grp == 0) | (grp == 2)
    stat, p = CX.logrank(t[is_hold][lo_hi], e[is_hold][lo_hi], (grp[lo_hi] == 2).astype(int))
    res["risk_groups"] = {"arm": best, "cutpoints_oof": [round(float(c), 4) for c in cut],
                          "km": km_curves, "logrank_low_vs_high": {"chi2": round(stat, 3), "p": p}}
    print("\n  risk tertiles on the hold-out (%s): low vs high log-rank chi2 %.2f, p = %.3g"
          % (best, stat, p))

    # ---------------- persist ----------------
    # BOTH deployable arms are kept, because an uploaded patient often has no age and no ELN. `full`
    # needs them; `molecular` does not but is meaningfully weaker (CV 0.689 vs 0.726). At inference the
    # layer picks the best arm it can actually feed and says which one it used, rather than silently
    # imputing a median age and reporting the strong model's accuracy.
    bundle = {"version": SurvivalModel.VERSION, "feature_space": fs, "horizons": HORIZONS,
              "arm_blocks": ARM_BLOCKS, "card": res, "models": {}, "risk_ref": {},
              # the exact mutation columns the model was fitted on, so inference can place the mutation
              # caller's predictions in the right slots instead of guessing at the layout
              "mut_columns": getattr(build_blocks, "mut_columns", []),
              # age-spline knots from the TRAINING quartiles, so inference rebuilds the same basis
              "age_knots": getattr(build_blocks, "age_knots", None),
              # symbol -> ENSG, so a sample keyed by gene symbol (the atlas) can be aligned at inference
              "sym2ens": __import__("train_drug_model").sym2ens_map(genes)}
    for arm in ("deployed", "full", "molecular", "clin", "age_eln"):
        m = fitted.get(arm)
        if m is None:
            continue
        bundle["models"][arm] = m
        Ball = {b: B[b] for b in ARM_BLOCKS[arm]}
        r_all = m.risk(Ball) if isinstance(m, SurvivalModel) else m.risk(B[ARM_BLOCKS[arm][0]])
        bundle["risk_ref"][arm] = np.sort(r_all[tr])
        if not isinstance(m, SurvivalModel):          # single-block arms need their own baseline hazard
            m._final = CX.CoxPH(alpha=1e-6).fit(r_all[tr].reshape(-1, 1), t[tr], e[tr])
    # --- treatment-stratified models ----------------------------------------------------------------
    # Measured, then deployed: a model fitted on the POOLED cohort and applied to non-intensively-treated
    # patients scores C-index 0.554, which is why this platform documented itself as "close to useless"
    # in that group. Fitting WITHIN the stratum gives 0.681. The failure was pooling, not the biology --
    # and the guideline fails there too (ELN 2022 scores 0.462, at or below chance), independently
    # reproducing Pollyea/Dohner (Blood 2024). See deliverables/ELN2022_RISK_BENCHMARK.md.
    #
    # Both strata are stored. Inference uses one only when the caller states the induction type; with no
    # treatment information the pooled model is still the right default, because guessing the stratum
    # would be worse than not stratifying.
    ty_all = cl["typeInductionTx"].astype(str)
    strata_masks = {"intensive": ty_all.str.contains("Standard Chemo", na=False).values}
    strata_masks["non_intensive"] = ~strata_masks["intensive"]
    bundle["strata"] = {}
    for sname, smask in strata_masks.items():
        sidx = np.where(smask & tr)[0]
        if len(sidx) < 40 or e[sidx].sum() < 15:
            print("  stratum %-14s SKIPPED (n=%d, %d events -- too few to fit)"
                  % (sname, int(smask.sum()), int(e[smask].sum())))
            continue
        try:
            _, Bs = build_blocks(cl, X_lin, mut_all, rows, sorted(rows[sidx]), genes,
                                 a.n_pc, a.var_genes, train_idx=sidx)
            ms = SurvivalModel().fit({b: Bs[b][sidx] for b in ARM_BLOCKS["deployed"]},
                                     t[sidx], e[sidx], g[sidx])
            rs = ms.risk({b: Bs[b] for b in ARM_BLOCKS["deployed"]})
            bundle["strata"][sname] = {"model": ms, "risk_ref": np.sort(rs[sidx]),
                                       "n_train": int(len(sidx)), "events_train": int(e[sidx].sum())}
            print("  stratum %-14s fitted on %d patients, %d events"
                  % (sname, len(sidx), int(e[sidx].sum())))
        except Exception as ex:
            print("  stratum %-14s FAILED: %s" % (sname, str(ex)[:90]))
    bundle["strata_note"] = ("stratum-specific deployed-arm models. Cross-validated C-index within "
                             "stratum: intensive 0.729, non-intensive 0.681, against 0.554 for the "
                             "pooled model applied to the non-intensive group.")

    # --- cohort-matched references for single-cell input -------------------------------------------
    # A single-cell bulk-equivalent z-scored against BeatAML lands far outside that distribution, which
    # drives the Cox linear predictor to an extreme and produces absurd curves (the first smoke test
    # told an atlas patient they had two weeks to live). Same fix as the drug layer: give the model an
    # `sc` expression reference and the sc cohort's own risk distribution, so a single-cell patient is
    # ranked among single-cell samples and only THEN mapped onto the BeatAML risk scale, whose baseline
    # hazard is the only one tied to actual observed survival.
    sc_X = d["sc_X"].astype(float)
    Zsc = fs.add_reference("sc", sc_X)
    Psc = fs.pca.transform(Zsc[:, fs.sel])
    Ssc, _ = fs._state_block(Zsc)
    ncol = len(getattr(build_blocks, "mut_columns", []) or [1])
    Msc = np.zeros((sc_X.shape[0], ncol)); Msc[:, -1] = 1.0        # no genotype panel for sc uploads
    sc_blocks = {"rna": Psc, "state": Ssc, "mut": Msc}
    for arm in ("molecular",):
        m = fitted.get(arm)
        if m is None:
            continue
        bundle["risk_ref_sc"] = {arm: np.sort(m.risk({b: sc_blocks[b] for b in ARM_BLOCKS[arm]}))}
    print("  sc reference: %d single-cell samples; median sc risk %.2f vs BeatAML %.2f"
          % (sc_X.shape[0], float(np.median(bundle["risk_ref_sc"]["molecular"])),
             float(np.median(bundle["risk_ref"]["molecular"]))))

    bundle["cohort_km"] = dict(zip(("times", "survival"),
                                   [[round(float(x), 3) for x in v] for v in CX.km(t[tr], e[tr])]))
    with open(OUTP, "wb") as f:
        pickle.dump(bundle, f)
    json.dump(res, open(CARD, "w"), indent=1)
    print("\nwrote %s + %s  (%.0fs)" % (OUTP, CARD, time.time() - t0))


if __name__ == "__main__":
    main()
