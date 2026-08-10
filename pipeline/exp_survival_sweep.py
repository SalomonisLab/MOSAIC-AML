#!/usr/bin/env python3
"""Roadmap B and C, run as one sweep. Every arm is patient-grouped CV on the same split as the deployed
model, so the numbers are directly comparable to C-index 0.726 (CV) / 0.752 (hold-out).

  B1  LSC17 stemness signature as benchmark and feature
  B3  variant detail currently discarded: FLT3-ITD allelic ratio, mutation VAF, karyotype features
  B4  ELN 2022 (derived earlier) vs ELN 2017
  B5  age as a spline rather than linear
  B6  elastic-net / univariate screening on RNA instead of 60 unscreened PCs
  B7  random survival forest and gradient-boosted survival
  L1  treatment-adjusted (time-dependent transplant, handled by landmark)
  L2  treatment-stratified (intensive vs not; transplant vs not)
  C8  proportional-hazards assumption test  -> decides whether C3/C4 are optional
  C3  accelerated failure time (log-normal / Weibull)
  C4  discrete-time multi-horizon classifiers
  C2  restricted mean survival time instead of median
  C1  conformal prediction intervals for survival time

  python exp_survival_sweep.py -> deliverables/exp_survival_sweep.json
"""
import os, sys, json, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score

from amlmm.survival import data as SD, coxph as CX
from amlmm.survival.model import SurvivalModel
from train_survival_model import build_blocks, fit_arm, ARM_BLOCKS, HORIZONS, BUNDLE, boot_delta_c

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_survival_sweep.json")
CLIN = os.path.join(ROOT, "data", "external", "beataml", "beataml_wv1to4_clinical.xlsx")
MUTS = os.path.join(ROOT, "data", "external", "beataml", "mutations.txt")
ELN22 = os.path.join(ROOT, "labels", "eln2022_beataml.tsv")

# LSC17 (Ng et al. 2016) — the published 17-gene AML stemness score, with its fitted weights.
LSC17 = {"DNMT3B": 0.0874, "ZBTB46": -0.0347, "NYNRIN": 0.00865, "ARHGAP22": -0.0138,
         "LAPTM4B": 0.00582, "MMRN1": 0.0258, "DPYSL3": 0.0284, "KIAA0125": 0.0196,
         "CDK6": -0.0704, "CPXM1": -0.0258, "SOCS2": 0.0271, "SMIM24": -0.0226,
         "EMP1": 0.0146, "NGFRAP1": 0.0465, "CD34": 0.0338, "AKR1C3": -0.0402, "GPR56": 0.0501}


def lsc17_score(X_lin, genes, sym2ens):
    """LSC17 on the shared ENSG space; returns (score, n_genes_found)."""
    gpos = {g: i for i, g in enumerate(genes)}
    L = np.log2(np.clip(X_lin, 0, None) + 1.0)
    Z = (L - L.mean(0)) / np.where(L.std(0) == 0, 1, L.std(0))
    s = np.zeros(X_lin.shape[0]); n = 0
    for sym, w in LSC17.items():
        e = sym2ens.get(sym)
        if e in gpos:
            s += w * Z[:, gpos[e]]; n += 1
    return s.reshape(-1, 1), n


def extra_clinical(cl):
    """B3 + B4: the prognostic detail currently thrown away."""
    cols, names = [], []
    ar = pd.to_numeric(cl.get("allelic_ratio"), errors="coerce")     # FLT3-ITD allelic ratio
    cols += [ar.fillna(0.0).values, ar.notna().astype(float).values,
             (ar.fillna(0) > 0.5).astype(float).values]
    names += ["flt3_itd_ar", "flt3_itd_ar_present", "flt3_itd_ar_high"]
    kar = cl.get("karyotype").astype(str).str.lower()
    ncl = kar.str.count(",")                                          # crude complexity proxy
    cols += [ncl.fillna(0).clip(0, 12).values,
             kar.str.contains("complex").astype(float).values,
             (ncl >= 5).astype(float).values,
             kar.str.contains(r"\-7|del\(7").astype(float).values,
             kar.str.contains(r"\-5|del\(5").astype(float).values,
             kar.str.contains(r"\+8").astype(float).values,
             kar.str.contains("46,x[xy]\\[").astype(float).values]
    names += ["kar_ncommas", "kar_complex_text", "kar_ge5_abn", "kar_del7", "kar_del5",
              "kar_tri8", "kar_normal"]
    fus = cl.get("consensusAMLFusions").astype(str)
    for f in ("RUNX1-RUNX1T1", "CBFB-MYH11", "PML-RARA", "KMT2A", "BCR-ABL1"):
        cols.append(fus.str.contains(f, case=False, na=False).astype(float).values)
        names.append("fus_" + f)
    e22 = {}
    if os.path.exists(ELN22):
        t = pd.read_csv(ELN22, sep="\t")
        t["specimen"] = t["dbgap_rnaseq_sample"].astype(str)
        e22 = dict(zip(t["specimen"], t["ELN2022"].astype(str)))
    v = cl["specimen"].map(lambda s: e22.get(s, ""))
    for lv in ("Favorable", "Intermediate", "Adverse"):
        cols.append(v.eq(lv).astype(float).values); names.append("eln2022_" + lv)
    cols.append(v.isin(["Favorable", "Intermediate", "Adverse"]).astype(float).values)
    names.append("eln2022_known")
    return np.vstack(cols).T.astype(float), names


def vaf_block(cl, muts_path):
    """Mutation VAF per driver gene — presence/absence discards clonal burden, which is prognostic."""
    GEN = ["TP53", "FLT3", "NPM1", "DNMT3A", "NRAS", "KRAS", "IDH1", "IDH2", "RUNX1", "ASXL1",
           "TET2", "SRSF2", "U2AF1", "STAG2", "CEBPA", "WT1", "PTPN11", "KIT", "EZH2", "BCOR"]
    m = pd.read_csv(muts_path, sep="\t", low_memory=False,
                    usecols=["dbgap_sample_id", "symbol", "t_vaf"])
    m["t_vaf"] = pd.to_numeric(m["t_vaf"], errors="coerce")
    top = m[m["symbol"].isin(GEN)].groupby(["dbgap_sample_id", "symbol"])["t_vaf"].max().unstack()
    dna = cl["dbgap_dnaseq_sample"].astype(str)
    T = top.reindex(dna.values)
    X = np.nan_to_num(T.reindex(columns=GEN).values.astype(float), nan=0.0)
    burden = np.nan_to_num(T.reindex(columns=GEN).notna().sum(1).values.astype(float))
    return np.hstack([X, burden.reshape(-1, 1)]), ["vaf_" + g for g in GEN] + ["n_drivers"]


def treatment_block(cl):
    """L1: treatment covariates. Transplant is deliberately NOT included as a baseline covariate —
    see the landmark analysis instead; coding it at baseline creates immortal-time bias."""
    cols, names = [], []
    ty = cl["typeInductionTx"].astype(str)
    cols += [ty.str.contains("Standard Chemo").astype(float).values,
             ty.str.contains("Targeted").astype(float).values,
             ty.str.contains("Palliative").astype(float).values]
    names += ["tx_intensive", "tx_targeted", "tx_palliative"]
    reg = cl["cumulativeTreatmentRegimens"].astype(str).str.lower()
    for tok in ("cytarabine", "azacitidine", "decitabine", "sorafenib", "midostaurin", "venetoclax"):
        cols.append(reg.str.contains(tok).astype(float).values); names.append("got_" + tok)
    cols.append(pd.to_numeric(cl.get("cumulativeTreatmentRegimenCount"), errors="coerce")
                .fillna(1).clip(0, 10).values)
    names.append("n_regimens")
    return np.vstack(cols).T.astype(float), names


# ------------------------------------------------------------------- C ----
def rmst(t, e, horizon):
    """Restricted mean survival time — defined even when the median is never reached."""
    ts, ss = CX.km(t, e)
    grid = np.concatenate([[0.0], ts[ts <= horizon], [horizon]])
    sv = np.concatenate([[1.0], ss[ts <= horizon], [np.interp(horizon, ts, ss, left=1.0,
                        right=ss[-1] if len(ss) else 1.0)]])
    return float(np.sum(np.diff(grid) * sv[:-1]))


def ph_test(X, t, e, beta):
    """Schoenfeld-residual correlation with time; a strong correlation means proportional hazards fails."""
    from scipy.stats import spearmanr
    o = np.lexsort((e == 0, t)); X, t, e = X[o], t[o], e[o]
    eta = X @ beta; w = np.exp(eta - eta.max())
    csum = np.cumsum(w[::-1])[::-1]
    wx = np.cumsum((w[:, None] * X)[::-1], axis=0)[::-1] / csum[:, None]
    ev = e.astype(bool)
    r = X[ev] - wx[ev]
    out = []
    for j in range(X.shape[1]):
        s = spearmanr(t[ev], r[:, j])
        out.append(float(s.pvalue))
    return {"n_covariates": X.shape[1], "min_p": round(float(np.min(out)), 5),
            "frac_p_lt_0.05": round(float(np.mean(np.asarray(out) < 0.05)), 3)}


def main():
    t0 = time.time()
    d = np.load(BUNDLE, allow_pickle=True)
    genes = [str(x) for x in d["genes"]]; ba = [str(s) for s in d["ba_samples"]]
    X_lin = d["ba_X"].astype(np.float64)
    mut_all = pd.DataFrame(d["ba_L"].astype(float), index=ba, columns=[str(c) for c in d["drivers"]])
    row_of = {s: i for i, s in enumerate(ba)}
    import train_drug_model as TD
    sym2ens = TD.sym2ens_map(genes)

    cl = SD.load_cohort(specimens=ba)
    rows = np.array([row_of[s] for s in cl["specimen"]])
    t = cl["time_years"].values.astype(float); e = cl["event"].values.astype(int)
    g = cl["subject"].values
    subs = np.array(sorted(set(g))); rng = np.random.RandomState(0)
    hold = set(subs[rng.permutation(len(subs))[:int(round(0.20 * len(subs)))]])
    is_hold = np.array([s in hold for s in g]); tr = ~is_hold; idx_tr = np.where(tr)[0]
    print("cohort %d patients, %d deaths | train %d, hold-out %d" % (len(t), e.sum(), tr.sum(), is_hold.sum()))

    # extra blocks
    L17, n17 = lsc17_score(X_lin[rows], genes, sym2ens)
    XC, ncl = extra_clinical(cl)
    XV, nvf = vaf_block(cl, MUTS)
    XT, ntx = treatment_block(cl)
    print("  LSC17 genes found %d/17 | extra clinical %d cols | VAF %d cols | treatment %d cols"
          % (n17, XC.shape[1], XV.shape[1], XT.shape[1]))

    ARMS = ["age_eln", "full", "lsc17", "full_lsc17", "full_variant", "full_treat", "full_all"]
    oof = {a: np.full(tr.sum(), np.nan) for a in ARMS}
    oof_extra = {a: np.full(tr.sum(), np.nan) for a in ("full_spline", "rna_screened", "gbm")}

    gk = GroupKFold(n_splits=5)
    for k, (i_in, i_out) in enumerate(gk.split(np.zeros(tr.sum()), groups=g[tr])):
        gi, go = idx_tr[i_in], idx_tr[i_out]
        _, B = build_blocks(cl, X_lin, mut_all, rows, sorted(rows[gi]), genes, 60, 4000)
        B["lsc17"] = L17; B["xtra"] = XC; B["vaf"] = XV; B["tx"] = XT
        # B5: age spline (natural cubic-ish via piecewise powers)
        age = B["age_eln"][:, [0]]
        kn = np.quantile(age[gi], [0.25, 0.5, 0.75])
        B["age_spline"] = np.hstack([age] + [np.clip(age - q, 0, None) ** 3 for q in kn])
        combos = {
            "age_eln": ["age_eln"], "full": ARM_BLOCKS["full"],
            "lsc17": ["lsc17"],
            "full_lsc17": ARM_BLOCKS["full"] + ["lsc17"],
            "full_variant": ARM_BLOCKS["full"] + ["lsc17", "xtra", "vaf"],
            "full_treat": ARM_BLOCKS["full"] + ["tx"],
            "full_all": ARM_BLOCKS["full"] + ["lsc17", "xtra", "vaf", "tx"],
            "full_spline": ["rna", "state", "mut", "clin", "age_spline"],
        }
        for name, blocks in combos.items():
            key = name if name in oof else name
            store = oof if name in oof else oof_extra
            try:
                if len(blocks) == 1:
                    m = CX.CoxPH(alpha=0.10).fit(B[blocks[0]][gi], t[gi], e[gi])
                    store[key][i_out] = m.risk(B[blocks[0]][go])
                else:
                    sm = SurvivalModel().fit({b: B[b][gi] for b in blocks}, t[gi], e[gi], g[gi])
                    store[key][i_out] = sm.risk({b: B[b][go] for b in blocks})
            except Exception as ex:
                print("   fold %d %s failed: %s" % (k + 1, name, ex))
        # B6: univariate-screened RNA
        try:
            P = B["rna"]
            sc = np.array([abs(CX.c_index(t[gi], e[gi], P[gi][:, j]) - 0.5) for j in range(P.shape[1])])
            sel = np.argsort(sc)[::-1][:15]
            m = CX.CoxPH(alpha=0.10).fit(P[gi][:, sel], t[gi], e[gi])
            oof_extra["rna_screened"][i_out] = m.risk(P[go][:, sel])
        except Exception as ex:
            print("   screening failed:", ex)
        # B7: gradient-boosted survival, approximated by GBM on the martingale residual
        try:
            from sklearn.ensemble import GradientBoostingRegressor
            base = CX.CoxPH(alpha=0.10).fit(B["clin"][gi], t[gi], e[gi])
            H0 = np.interp(t[gi], base.baseline_t, base.baseline_H, left=0,
                           right=base.baseline_H[-1] if len(base.baseline_H) else 0)
            mart = e[gi] - H0 * np.exp(base.risk(B["clin"][gi]))
            F = np.hstack([B["rna"], B["state"], B["clin"], B["lsc17"], B["xtra"], B["vaf"]])
            gb = GradientBoostingRegressor(n_estimators=250, max_depth=2, learning_rate=0.05,
                                           subsample=0.8, random_state=0).fit(F[gi], mart)
            # a LARGE positive martingale residual = died when the model expected survival = higher
            # risk, so the predicted residual IS the risk (the earlier minus sign inverted it)
            oof_extra["gbm"][i_out] = gb.predict(F[go])
        except Exception as ex:
            print("   gbm failed:", ex)
        print("   fold %d/5 (%.0fs)" % (k + 1, time.time() - t0))

    t_tr, e_tr = t[tr], e[tr]
    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "n": int(tr.sum()), "events": int(e_tr.sum()),
           "lsc17_genes_found": n17, "arms": {}}
    allp = dict(oof); allp.update(oof_extra)
    base = allp["age_eln"]
    print("\n== B · survival discrimination (patient-grouped CV, n=%d, %d deaths) ==" % (tr.sum(), e_tr.sum()))
    print("  %-14s %8s %10s" % ("arm", "C-index", "vs age+ELN"))
    for a, r in allp.items():
        m = np.isfinite(r)
        if m.sum() < 50:
            print("  %-14s   (not scored)" % a); continue
        c = CX.c_index(t_tr[m], e_tr[m], r[m])
        dlt = boot_delta_c(t_tr, e_tr, r, base, B=400) if a != "age_eln" else None
        res["arms"][a] = {"c_index": round(c, 4), "n_scored": int(m.sum()),
                          "delta_vs_age_eln": dlt}
        print("  %-14s %8.3f %10s" % (a, c, "" if dlt is None else
              "%+.4f [%+.3f,%+.3f]" % (dlt["delta_c"], dlt["ci95"][0], dlt["ci95"][1])))

    # ---- L2 treatment-stratified ----
    print("\n== L2 · treatment-stratified ==")
    intensive = XT[:, 0][tr].astype(bool)
    strat = {}
    for nm, msk in (("intensive induction", intensive), ("non-intensive", ~intensive)):
        r = allp["full"]; m = msk & np.isfinite(r)
        if m.sum() >= 40 and e_tr[m].sum() >= 15:
            strat[nm] = {"n": int(m.sum()), "events": int(e_tr[m].sum()),
                         "c_full": round(CX.c_index(t_tr[m], e_tr[m], r[m]), 4),
                         "c_age_eln": round(CX.c_index(t_tr[m], e_tr[m], base[m]), 4)}
            print("  %-22s n=%-4d deaths=%-4d C full %.3f vs age+ELN %.3f"
                  % (nm, strat[nm]["n"], strat[nm]["events"], strat[nm]["c_full"], strat[nm]["c_age_eln"]))
    res["L2_treatment_stratified"] = strat

    # ---- C8 proportional hazards ----
    Xph = np.column_stack([allp["full"][np.isfinite(allp["full"])]])
    mph = np.isfinite(allp["full"])
    cph = CX.CoxPH(alpha=1e-6).fit(Xph, t_tr[mph], e_tr[mph])
    res["C8_proportional_hazards"] = ph_test(Xph, t_tr[mph], e_tr[mph], cph.beta)
    print("\n== C8 · proportional-hazards test == Schoenfeld min p = %s"
          % res["C8_proportional_hazards"]["min_p"])

    # ---- C2 RMST + C3 AFT + C1 conformal, on the risk score ----
    r = allp["full"]; m = np.isfinite(r)
    cut = np.quantile(r[m], [1/3, 2/3]); grp = np.digitize(r[m], cut)
    res["C2_rmst"] = {}
    for gi_, nm in enumerate(("low", "intermediate", "high")):
        s = grp == gi_
        if s.sum() >= 10:
            res["C2_rmst"][nm] = {"n": int(s.sum()),
                                  "rmst_2y": round(rmst(t_tr[m][s], e_tr[m][s], 2.0), 3),
                                  "rmst_5y": round(rmst(t_tr[m][s], e_tr[m][s], 5.0), 3)}
    print("\n== C2 · restricted mean survival time (defined even when the median is not reached) ==")
    for k_, v in res["C2_rmst"].items():
        print("  %-13s n=%-4d RMST(2y) %.2f y | RMST(5y) %.2f y" % (k_, v["n"], v["rmst_2y"], v["rmst_5y"]))

    # C3 AFT (log-normal, censoring-aware) on the same features, compared by C-index on -predicted time
    try:
        from scipy.optimize import minimize
        from scipy.stats import norm
        Xa = np.column_stack([np.ones(m.sum()), r[m]])
        y = np.log(np.clip(t_tr[m], 1e-3, None)); dd = e_tr[m]
        def nll(p):
            b, ls = p[:-1], p[-1]; s = np.exp(ls); z = (y - Xa @ b) / s
            return -np.sum(dd * (norm.logpdf(z) - ls) + (1 - dd) * norm.logsf(z))
        opt = minimize(nll, np.r_[np.zeros(Xa.shape[1]), 0.0], method="Nelder-Mead",
                       options={"maxiter": 4000})
        pred_log_t = Xa @ opt.x[:-1]
        res["C3_aft_lognormal"] = {"c_index": round(CX.c_index(t_tr[m], e_tr[m], -pred_log_t), 4),
                                   "sigma": round(float(np.exp(opt.x[-1])), 3)}
        print("\n== C3 · log-normal AFT == C-index %.3f (Cox %.3f), sigma %.2f"
              % (res["C3_aft_lognormal"]["c_index"], res["arms"]["full"]["c_index"],
                 res["C3_aft_lognormal"]["sigma"]))
    except Exception as ex:
        print("  AFT failed:", ex)

    # C1 conformal interval on log-time, calibrated on held-out deaths
    try:
        died = (e_tr[m] == 1)
        resid = np.abs(y[died] - pred_log_t[died])
        for cov in (0.8, 0.9):
            q = float(np.quantile(resid, cov))
            lo = np.exp(pred_log_t - q); hi = np.exp(pred_log_t + q)
            emp = float(np.mean((t_tr[m][died] >= lo[died]) & (t_tr[m][died] <= hi[died])))
            res.setdefault("C1_conformal", {})["%.0f%%" % (100 * cov)] = {
                "target_coverage": cov, "empirical_coverage_on_deaths": round(emp, 3),
                "median_interval_width_years": round(float(np.median(hi - lo)), 2)}
        print("\n== C1 · conformal prediction intervals ==")
        for k_, v in res["C1_conformal"].items():
            print("  target %s -> empirical %.0f%%, median width %.1f y"
                  % (k_, 100 * v["empirical_coverage_on_deaths"], v["median_interval_width_years"]))
    except Exception as ex:
        print("  conformal failed:", ex)

    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
