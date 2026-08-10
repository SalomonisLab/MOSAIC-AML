"""The treatment-utility score, and why it is not just "rank by predicted sensitivity".

    S_ij = w1*P(sensitive) + w2*coverage_blast + w3*coverage_LSC + w4*mechanistic + w5*clinical
           - w6*uncertainty - w7*resistance - w8*infeasibility - w9*OOD

Every negative term exists because some real failure mode would otherwise rank first:

  uncertainty     a drug the model cannot predict at all (OOF AUROC ~0.5) must not outrank one it
                  predicts well, however extreme its point estimate
  resistance      a bulk-sensitive prediction with a resistant primitive subclone is a relapse, not a
                  response
  infeasibility   an intensive agent in an 80-year-old is not a recommendation
  OOD             a patient outside the training distribution gets a discounted score, not a confident
                  one, because the model has no basis for the extrapolation

Rankings are produced SEPARATELY per clinical tier. A research-only tool compound that scores 0.9 and
an approved AML therapy that scores 0.6 are not comparable options, and a single merged list invites
exactly the misreading the report exists to prevent. Nothing here is a treatment recommendation: the
output is a prioritisation for trial matching or laboratory validation.
"""
from __future__ import annotations
import numpy as np

from . import targets as TG

WEIGHTS = {"sensitivity": 0.30, "coverage_blast": 0.15, "coverage_lsc": 0.15,
           "mechanistic": 0.15, "clinical": 0.15,
           "uncertainty": 0.15, "resistance": 0.15, "infeasibility": 0.10, "ood": 0.20}


def uncertainty(prob, oof_auroc, n_train, spread=None, assay_reliability=None):
    """Four independent reasons to be unsure: can this drug be predicted at all, is this particular
    call decisive, how much data stands behind it — and does the underlying ASSAY even reproduce.

    The last one was added after measuring it: 19 of 79 inhibitors have a same-mechanism reliability
    below 0.15, and for those the model already sits at the achievable ceiling. A confident-looking
    recommendation for an inhibitor whose own measurement does not reproduce is the most misleading
    output this layer can produce, so low reliability now raises uncertainty directly."""
    u_model = 1.0 - max(0.0, min(1.0, ((oof_auroc or 0.5) - 0.5) * 2.0))
    u_point = 1.0 if prob is None else 1.0 - 2.0 * abs(float(prob) - 0.5)
    u_data = 1.0 / (1.0 + (n_train or 0) / 200.0)
    parts = [u_model, u_point, u_data]
    w = [0.45, 0.35, 0.20]
    if assay_reliability is not None and assay_reliability == assay_reliability:
        # reliability 0.45+ -> no penalty; 0.0 -> full penalty
        parts.append(float(min(1.0, max(0.0, 1.0 - assay_reliability / 0.45)))); w.append(0.40)
    if spread is not None:                      # across-cell-state spread, when Model B ran
        parts.append(min(1.0, float(spread) / 1.5)); w.append(0.25)
    return float(np.average(parts, weights=w[:len(parts)]))


def resistance_burden(mech_evidence, state_metrics=None, prob_cut=0.5):
    """Flagged bypass mechanisms, plus an explicit escape term: a resistant primitive subclone under a
    sensitive bulk prediction is the classic pattern behind ex-vivo-response-guided relapse."""
    flags = (mech_evidence or {}).get("resistance") or []
    meas = [r for r in flags if r.get("measurable")]
    frac = (sum(1 for r in meas if r.get("flagged")) / len(meas)) if meas else 0.0
    escape = 0.0
    if state_metrics:
        cov = state_metrics.get("coverage_LSC_like")
        pb = state_metrics.get("prob_sensitive_bulk")
        if cov is not None and pb is not None and pb >= prob_cut:
            escape = float(max(0.0, 1.0 - cov))        # bulk says yes, the primitive compartment does not
    return float(min(1.0, 0.6 * frac + 0.6 * escape)), {"flagged_fraction": round(frac, 3),
                                                        "subclone_escape": round(escape, 3),
                                                        "n_measurable": len(meas)}


def infeasibility(drug, clinical=None):
    """Deliberately conservative and explicit about what it does not know. With no age, performance
    status or organ function we return 0 and say 'not assessed' rather than inventing a penalty."""
    if not clinical:
        return 0.0, {"assessed": False, "reason": "no clinical covariates supplied"}
    pen, why = 0.0, []
    age = clinical.get("age")
    tier = TG.get(drug)["clinical_tier"]
    if age is not None and age >= 75 and tier == "research":
        pen += 0.3; why.append("research compound in a patient >=75")
    if age is not None and age >= 75:
        pen += 0.1; why.append("age >= 75")
    if clinical.get("prior_lines") and clinical["prior_lines"] >= 3 and tier == "research":
        pen += 0.2; why.append("heavily pre-treated, research-only compound")
    return float(min(1.0, pen)), {"assessed": True, "reasons": why}


def ood_penalty(distance, ref_distances):
    """Where this patient sits in the training cohort's own distance distribution. A patient at the
    cohort median gets ~0; one beyond everything the model was fitted on approaches 1."""
    if distance is None or ref_distances is None or not len(ref_distances):
        return 0.0
    q = float((np.asarray(ref_distances) <= float(distance)).mean())
    return float(max(0.0, (q - 0.5) / 0.5))          # no penalty until above the cohort median


def score(drug, prob_sensitive, oof_auroc, n_train, mech_evidence=None, state_metrics=None,
          clinical=None, ood_distance=None, ood_reference=None, weights=None,
          assay_reliability=None):
    w = dict(WEIGHTS); w.update(weights or {})
    sm = state_metrics or {}
    cov_b = sm.get("coverage_blast")
    cov_l = sm.get("coverage_LSC_like")
    mech = (mech_evidence or {}).get("mechanistic_score")
    clin_s = TG.clinical_score(drug)
    u = uncertainty(prob_sensitive, oof_auroc, n_train, sm.get("dispersion"), assay_reliability)
    r, r_detail = resistance_burden(mech_evidence, sm)
    f, f_detail = infeasibility(drug, clinical)
    o = ood_penalty(ood_distance, ood_reference)

    pos = (w["sensitivity"] * (0.0 if prob_sensitive is None else float(prob_sensitive))
           + w["coverage_blast"] * (0.0 if cov_b is None else cov_b)
           + w["coverage_lsc"] * (0.0 if cov_l is None else cov_l)
           + w["mechanistic"] * (0.0 if mech is None else mech)
           + w["clinical"] * clin_s)
    neg = (w["uncertainty"] * u + w["resistance"] * r + w["infeasibility"] * f + w["ood"] * o)
    # normalise the positive part by the weights that were actually evaluable, so a patient without
    # single-cell data is not silently penalised for missing coverage terms
    avail = (w["sensitivity"] * (prob_sensitive is not None) + w["coverage_blast"] * (cov_b is not None)
             + w["coverage_lsc"] * (cov_l is not None) + w["mechanistic"] * (mech is not None)
             + w["clinical"])
    pos = pos / avail if avail > 0 else 0.0
    return {
        "utility": round(float(pos - neg), 4),
        "components": {
            "sensitivity": None if prob_sensitive is None else round(float(prob_sensitive), 4),
            "coverage_blast": cov_b, "coverage_LSC_like": cov_l,
            "mechanistic": mech, "clinical_evidence": round(clin_s, 3),
            "uncertainty": round(u, 4), "assay_reliability": assay_reliability,
            "resistance": round(r, 4),
            "infeasibility": round(f, 4), "ood": round(o, 4),
        },
        "detail": {"resistance": r_detail, "infeasibility": f_detail,
                   "evaluable_positive_weight": round(float(avail), 3)},
        "weights": w,
    }


def rank_by_tier(scored, min_utility=None):
    """Four ordered lists, never one. `scored` = {drug: score(...) dict}."""
    out = {}
    for tier in TG.TIERS:
        rows = [{"inhibitor": d, **v} for d, v in scored.items()
                if TG.get(d)["clinical_tier"] == tier
                and (min_utility is None or v["utility"] >= min_utility)]
        rows.sort(key=lambda r: -r["utility"])
        out[tier] = {"label": TG.TIER_LABEL[tier], "n": len(rows), "ranked": rows}
    return out
