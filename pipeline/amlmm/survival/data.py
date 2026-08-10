"""The survival cohort: who can honestly be used to learn survival-from-diagnosis, and who cannot.

Three restrictions, each of which changes the answer:

1. **Known vital status.** 22 BeatAML specimens with expression have `vitalStatus == Unknown`. They have
   a follow-up time but no outcome, so they are neither events nor censorings and are dropped.
2. **Initial-diagnosis specimens only.** `overallSurvival` is measured from diagnosis. Using a relapse
   or residual-disease specimen's expression to predict survival-from-diagnosis mixes two different
   questions and quietly leaks prognosis (a patient with a relapse specimen has, by definition, already
   survived to relapse). That takes 649 usable specimens down to 446 — worth it.
3. **One specimen per patient.** Two patients contributed more than one initial-diagnosis specimen;
   the earliest is kept so that a patient appears exactly once.

Censoring is respected throughout: `event = 1` for a death, `0` for a patient last known alive, whose
`overallSurvival` is follow-up time rather than survival time.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
BA_DIR = os.path.join(ROOT, "data", "external", "beataml")
CLIN = "beataml_wv1to4_clinical.xlsx"

DAYS_PER_YEAR = 365.25


def load_cohort(specimens=None, ba_dir=BA_DIR, initial_only=True):
    """-> DataFrame with specimen, subject, time_years, event, and the clinical covariates."""
    cl = pd.read_excel(os.path.join(ba_dir, CLIN))
    cl["specimen"] = cl["dbgap_rnaseq_sample"].astype(str)
    cl["subject"] = cl["dbgap_subject_id"].astype(str)
    if specimens is not None:
        cl = cl[cl["specimen"].isin(set(map(str, specimens)))]
    n_all = len(cl)

    os_days = pd.to_numeric(cl["overallSurvival"], errors="coerce")
    keep = cl["vitalStatus"].isin(["Dead", "Alive"]) & os_days.notna() & (os_days >= 0)
    dropped_unknown = int((~cl["vitalStatus"].isin(["Dead", "Alive"])).sum())
    cl = cl[keep].copy()
    cl["time_years"] = os_days[keep] / DAYS_PER_YEAR
    cl["event"] = (cl["vitalStatus"] == "Dead").astype(int)

    n_before_stage = len(cl)
    if initial_only:
        cl = cl[cl["diseaseStageAtSpecimenCollection"].eq("Initial Diagnosis")].copy()
    # one specimen per patient: earliest collection
    cl["_order"] = pd.to_numeric(cl.get("timeOfSampleCollectionRelativeToInclusion"),
                                 errors="coerce").fillna(0)
    cl = cl.sort_values("_order").drop_duplicates("subject", keep="first").drop(columns="_order")

    cl.attrs["provenance"] = {
        "clinical_rows_with_expression": int(n_all),
        "dropped_unknown_vital_status": dropped_unknown,
        "after_outcome_filter": int(n_before_stage),
        "initial_diagnosis_only": bool(initial_only),
        "final_patients": int(len(cl)),
        "events": int(cl["event"].sum()),
        "censored": int((cl["event"] == 0).sum()),
        "median_followup_years_censored": round(float(cl.loc[cl["event"] == 0, "time_years"].median()), 2)
                                          if (cl["event"] == 0).any() else None,
    }
    return cl.reset_index(drop=True)


ELN_ORDER = {"Favorable": 0, "FavorableOrIntermediate": 0.5, "Intermediate": 1,
             "IntermediateOrAdverse": 1.5, "Adverse": 2}


def clinical_block(cl):
    """The baseline every molecular model has to beat: age and ELN 2017, plus routine labs.

    ELN enters both as an ordinal severity and as indicator columns, because the ordinal encoding
    assumes equal spacing between risk tiers and the indicators do not; the ridge decides.
    """
    cols, names = [], []
    age = pd.to_numeric(cl["ageAtDiagnosis"], errors="coerce")
    cols.append(age.fillna(age.median()).values); names.append("age")
    cols.append(age.isna().astype(float).values); names.append("age__missing")
    eln = cl["ELN2017"].map(ELN_ORDER)
    cols.append(eln.fillna(1.0).values); names.append("eln_ordinal")
    for lv in ("Favorable", "Intermediate", "Adverse"):
        cols.append(cl["ELN2017"].eq(lv).astype(float).values); names.append("eln_" + lv)
    cols.append(cl["consensus_sex"].eq("Male").astype(float).values); names.append("sex_male")
    for c, nm in ((("%.Blasts.in.BM"), "blasts_bm"), (("%.Blasts.in.PB"), "blasts_pb"),
                  ("wbcCount", "wbc")):
        v = pd.to_numeric(cl.get(c), errors="coerce")
        cols.append(v.fillna(v.median()).values); names.append(nm)
        cols.append(v.isna().astype(float).values); names.append(nm + "__missing")
    for c, nm in (("isRelapse", "is_relapse"), ("isDenovo", "is_denovo"),
                  ("isTransformed", "is_transformed"), ("priorMDS", "prior_mds")):
        s = cl.get(c)
        cols.append(pd.Series(s).astype(str).str.lower().isin(["true", "y", "yes", "1"])
                    .astype(float).values if s is not None else np.zeros(len(cl)))
        names.append(nm)
    X = np.vstack(cols).T.astype(float)
    return X, names


def age_eln_block(cl):
    """The narrow clinical baseline used for the head-to-head: age + ELN only."""
    age = pd.to_numeric(cl["ageAtDiagnosis"], errors="coerce")
    eln = cl["ELN2017"].map(ELN_ORDER)
    X = np.vstack([age.fillna(age.median()).values, eln.fillna(1.0).values,
                   cl["ELN2017"].eq("Adverse").astype(float).values,
                   cl["ELN2017"].eq("Favorable").astype(float).values]).T
    return X.astype(float), ["age", "eln_ordinal", "eln_Adverse", "eln_Favorable"]
