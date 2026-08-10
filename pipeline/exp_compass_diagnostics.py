#!/usr/bin/env python3
"""Roadmap A2 + A1 — the two experiments that have to run before any tuning.

A2. Is the per-drug AUROC partly a PATIENT main effect?
    Some specimens' cells simply die readily in culture (blast fraction, viability, handling). A model
    that learns "this specimen is globally sensitive" scores well on every drug while knowing nothing
    drug-specific. Decompose the response matrix into patient + drug + interaction, then re-score the
    deployed predictions against the INTERACTION ALONE. That is the honest drug-specific number.

A1. Does COMPASS predict PATIENTS, not just the assay?
    254 patients have an inhibitor screen, a clean complete-response/refractory label and survival.
    Ask whether predicted ex-vivo sensitivity to the drug they ACTUALLY RECEIVED tracks their clinical
    response and their overall survival. Out-of-fold predictions only.

  python exp_compass_diagnostics.py -> deliverables/exp_compass_diagnostics.json
"""
import os, sys, json, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from amlmm.drug import data as D, targets as TG
from amlmm.survival import coxph as CX
from eval_drug_model import load_all, oof_predictions, per_drug_metrics, per_patient_metrics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_compass_diagnostics.json")
CLIN = os.path.join(ROOT, "data", "external", "beataml", "beataml_wv1to4_clinical.xlsx")


# --------------------------------------------------------------------- A2 ----
def two_way_decompose(df):
    """M_ij = mu + patient_i + drug_j + interaction_ij, by alternating median centring (robust to the
    skewed AUC distribution). Returns the components and how much variance each explains."""
    M = df.pivot_table(index="specimen", columns="inhibitor", values="auc", aggfunc="mean")
    obs = M.notna()
    mu = float(np.nanmedian(M.values))
    R = M - mu
    pat = pd.Series(0.0, index=M.index); drg = pd.Series(0.0, index=M.columns)
    for _ in range(50):
        dp = R.median(axis=1).fillna(0.0)                 # patient effect
        R = R.sub(dp, axis=0); pat += dp
        dd = R.median(axis=0).fillna(0.0)                 # drug effect
        R = R.sub(dd, axis=1); drg += dd
        if float(dp.abs().max()) < 1e-6 and float(dd.abs().max()) < 1e-6:
            break
    tot = np.nanvar(M.values[obs.values])
    v_pat = np.nanvar(np.add.outer(pat.values, np.zeros(len(drg)))[obs.values])
    v_drg = np.nanvar(np.add.outer(np.zeros(len(pat)), drg.values)[obs.values])
    v_int = np.nanvar(R.values[obs.values])
    return {"matrix": M, "resid": R, "patient": pat, "drug": drg, "mu": mu,
            "variance": {"total": float(tot),
                         "patient_effect": round(float(v_pat / tot), 4),
                         "drug_effect": round(float(v_drg / tot), 4),
                         "interaction": round(float(v_int / tot), 4)}}


def a2(A, df):
    """Re-score the deployed OOF predictions against the double-centred (interaction) target."""
    dec = two_way_decompose(A["long"])
    R = dec["resid"]
    key = list(zip(df["specimen"], df["inhibitor"]))
    inter = np.array([R.at[s, d] if (s in R.index and d in R.columns) else np.nan for s, d in key])
    d2 = df.copy()
    d2["sens_interaction"] = -inter                       # sign: higher = more sensitive, as elsewhere

    rows = []
    for drug, g in d2.groupby("inhibitor"):
        m = g["pred"].notna() & g["sens"].notna() & g["sens_interaction"].notna()
        if int(m.sum()) < 30:
            continue
        y_cur = g.loc[m, "sens"].values                    # within-drug z  = patient effect + interaction
        y_int = g.loc[m, "sens_interaction"].values        # interaction only
        p = g.loc[m, "pred"].values
        # tail AUROC on the interaction target, tails defined the same way (bottom/top 20%)
        lo, hi = np.quantile(y_int, [0.2, 0.8])
        yy = np.where(y_int >= hi, 1, np.where(y_int <= lo, 0, -1))
        k = yy >= 0
        rows.append({"inhibitor": drug, "n": int(m.sum()),
                     "spearman_current_target": round(float(spearmanr(y_cur, p).statistic), 4),
                     "spearman_interaction_only": round(float(spearmanr(y_int, p).statistic), 4),
                     "auroc_current_target": None,
                     "auroc_interaction_only": (round(float(roc_auc_score(yy[k], p[k])), 4)
                                                if k.sum() >= 20 and len(set(yy[k])) == 2 else None)})
    t = pd.DataFrame(rows)
    # how well does the model predict the PATIENT main effect itself?
    pat = dec["patient"]
    per_spec = d2.groupby("specimen")["pred"].mean()
    common = [s for s in per_spec.index if s in pat.index]
    r_pat = spearmanr(per_spec[common].values, -pat[common].values)
    return {
        "variance_decomposition": dec["variance"],
        "n_drugs_scored": int(len(t)),
        "mean_spearman_current_target": round(float(t["spearman_current_target"].mean()), 4),
        "mean_spearman_interaction_only": round(float(t["spearman_interaction_only"].mean()), 4),
        "mean_auroc_interaction_only": round(float(t["auroc_interaction_only"].mean()), 4),
        "model_vs_patient_main_effect_spearman": round(float(r_pat.statistic), 4),
        "interpretation": None,
        "table": t.sort_values("spearman_interaction_only", ascending=False).head(20).to_dict("records"),
    }, dec


# --------------------------------------------------------------------- A1 ----
DRUG_FOR_REGIMEN = {"cytarabine": "Cytarabine", "azacitidine": "Azacytidine",
                    "decitabine": "Azacytidine",            # nearest modelled HMA
                    "sorafenib": "Sorafenib", "midostaurin": "Midostaurin",
                    "venetoclax": "Venetoclax"}


def a1(df):
    """Does OOF-predicted ex-vivo sensitivity to the drug the patient ACTUALLY GOT track their outcome?"""
    cl = pd.read_excel(CLIN)
    cl["specimen"] = cl["dbgap_rnaseq_sample"].astype(str)
    cl = cl[cl["diseaseStageAtSpecimenCollection"].eq("Initial Diagnosis")].drop_duplicates("specimen")
    r = cl["responseToInductionTx"].astype(str)
    cl["cr"] = np.where(r.str.startswith("Complete Response"), 1,
                        np.where(r.eq("Refractory"), 0, np.nan))
    cl["os_years"] = pd.to_numeric(cl["overallSurvival"], errors="coerce") / 365.25
    cl["dead"] = (cl["vitalStatus"] == "Dead").astype(int)
    cl["age"] = pd.to_numeric(cl["ageAtDiagnosis"], errors="coerce")
    eln_ord = {"Favorable": 0, "Intermediate": 1, "Adverse": 2}
    cl["eln"] = cl["ELN2017"].map(eln_ord)
    reg = cl["cumulativeTreatmentRegimens"].astype(str).str.lower()

    pred = df.pivot_table(index="specimen", columns="inhibitor", values="pred", aggfunc="mean")
    out = {}
    for token, drug in DRUG_FOR_REGIMEN.items():
        if drug not in pred.columns:
            continue
        got = cl[reg.str.contains(token, na=False)].copy()
        got["p"] = got["specimen"].map(pred[drug])
        g = got[got["p"].notna()]
        e = {"drug_model": drug, "regimen_token": token, "n_received_and_scored": int(len(g))}
        # (a) clinical response
        gg = g[g["cr"].notna()]
        if len(gg) >= 25 and gg["cr"].nunique() == 2:
            e["response"] = {
                "n": int(len(gg)), "n_CR": int(gg["cr"].sum()),
                "auroc_predicted_sensitivity_vs_CR": round(float(roc_auc_score(gg["cr"], gg["p"])), 4),
                "spearman": round(float(spearmanr(gg["cr"], gg["p"]).statistic), 4),
                "p_value": round(float(spearmanr(gg["cr"], gg["p"]).pvalue), 4)}
        # (b) overall survival, unadjusted and adjusted for age + ELN
        gs = g[g["os_years"].notna()]
        if len(gs) >= 30 and gs["dead"].sum() >= 15:
            t_, e_, x = gs["os_years"].values, gs["dead"].values, gs["p"].values
            e["survival"] = {"n": int(len(gs)), "deaths": int(e_.sum()),
                             "c_index_predicted_sensitivity": round(CX.c_index(t_, e_, -x), 4)}
            adj = gs.dropna(subset=["age", "eln"])
            if len(adj) >= 30 and adj["dead"].sum() >= 15:
                Xa = np.column_stack([adj["age"].values, adj["eln"].values])
                Xb = np.column_stack([adj["age"].values, adj["eln"].values, adj["p"].values])
                ta, ea = adj["os_years"].values, adj["dead"].values
                m1 = CX.CoxPH(alpha=0.01).fit(Xa, ta, ea)
                m2 = CX.CoxPH(alpha=0.01).fit(Xb, ta, ea)
                e["survival"].update({
                    "n_adjusted": int(len(adj)),
                    "c_age_eln": round(CX.c_index(ta, ea, m1.risk(Xa)), 4),
                    "c_age_eln_plus_predicted": round(CX.c_index(ta, ea, m2.risk(Xb)), 4),
                    "beta_predicted_sensitivity": round(float(m2.beta[-1]), 4),
                    "note": "negative beta = predicted-sensitive patients died LESS often (expected direction)"})
        out[drug + " (" + token + ")"] = e
    return out


def main():
    t0 = time.time()
    A = load_all()
    print("scoring the cohort out-of-fold (this is the same 5-fold donor-grouped pass as the eval)…")
    df = oof_predictions(A, folds=5)
    print("  done (%.0fs)" % (time.time() - t0))

    res = {"generated": time.strftime("%Y-%m-%d %H:%M")}
    print("\n=== A2 · is the headline a patient main effect? ===")
    r2, dec = a2(A, df)
    v = r2["variance_decomposition"]
    print("  variance of the AUC matrix: patient %.1f%% | drug %.1f%% | interaction %.1f%%"
          % (100 * v["patient_effect"], 100 * v["drug_effect"], 100 * v["interaction"]))
    print("  mean per-drug Spearman  current target %.3f  ->  interaction only %.3f"
          % (r2["mean_spearman_current_target"], r2["mean_spearman_interaction_only"]))
    print("  mean per-drug AUROC on the interaction-only target: %.3f" % r2["mean_auroc_interaction_only"])
    print("  model's per-patient mean prediction vs the patient main effect: Spearman %.3f"
          % r2["model_vs_patient_main_effect_spearman"])
    keep = r2["mean_spearman_interaction_only"] / max(1e-9, r2["mean_spearman_current_target"])
    r2["interpretation"] = ("%.0f%% of the current per-drug Spearman survives removing the patient main "
                            "effect" % (100 * keep))
    print("  -> " + r2["interpretation"])
    res["A2_patient_vs_drug_effect"] = r2

    print("\n=== A1 · does it predict PATIENTS, not just the assay? ===")
    r1 = a1(df)
    for k, e in r1.items():
        line = "  %-28s n=%-4d" % (k, e["n_received_and_scored"])
        if "response" in e:
            line += " | CR AUROC %.3f (n=%d, %d CR, p=%.3f)" % (
                e["response"]["auroc_predicted_sensitivity_vs_CR"], e["response"]["n"],
                e["response"]["n_CR"], e["response"]["p_value"])
        if "survival" in e and "c_age_eln" in e["survival"]:
            s = e["survival"]
            line += " | OS C %.3f -> %.3f with model" % (s["c_age_eln"], s["c_age_eln_plus_predicted"])
        print(line)
    res["A1_clinical_outcome_validation"] = r1

    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
