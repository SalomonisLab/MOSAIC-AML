"""The eight drug-reasoning expert agents.

Same contract as the rest of the MOSAIC panel: each agent is a DETERMINISTIC evidence gatherer with a
declared `grounding` and `independence`, and the LLM (if one is attached) narrates within conclusions
it cannot change. All eight are `therapeutic`-domain and therefore **non-voting** with respect to the
anchored subtype call -- a drug recommendation must never be able to reach back and alter the genetic
anchor. They produce a separate output: a tiered, evidence-carrying prioritisation.

  drug_response        predicted response, calibrated probability, percentile, uncertainty, and the
                       nearest BeatAML specimens with their MEASURED AUCs
  cell_state_coverage  which cell states respond, which escape, and how much of the blast and LSC-like
                       compartment each drug actually covers
  molecular_mechanism  Model C: target, pathway readout, dependency, genetic activation, bypass routes
  pharmacology         whether the apparent sensitivity lives at concentrations plausibly reachable in
                       a patient, and whether it is on-target
  clinical_evidence    tier separation and eligibility, so a tool compound cannot masquerade as therapy
  combination          coverage-complementarity hypotheses only -- BeatAML measured single agents, so
                       additive single-agent scores are NOT evidence of synergy and are refused
  skeptic              nine specific challenges, each of which can demote a recommendation
  reporting            assembles the structured recommendation; refuses to emit a narrative
"""
from __future__ import annotations
import numpy as np

from . import targets as TG, utility as U, statemodel as SM

DOMAIN = "therapeutic"          # non-voting: cannot influence the anchored subtype call


def _res(name, grounding, independence, evidence, opinion=None):
    return {"name": name, "domain": DOMAIN, "grounding": grounding, "independence": independence,
            "status": "ok", "evidence": evidence, "opinion": opinion or {}, "votes": False}


# ------------------------------------------------------------ 1. response ----
def drug_response_agent(mod, preds, neighbours=None, top=15):
    """Predicted response with its own honest reliability, plus real measured AUCs from the most
    similar training specimens -- the part of the evidence a reader can check without the model."""
    oof = getattr(mod, "oof_metrics", {}) or {}
    rows = []
    for drug, p in preds.items():
        m = oof.get(drug, {})
        row = {"inhibitor": drug, "sens": p.get("sens"),
               "prob_sensitive": p.get("prob_sensitive"), "percentile": p.get("percentile"),
               "n_train": p.get("n_train"),
               "model_oof_auroc": m.get("auroc"), "model_oof_spearman": m.get("spearman"),
               "uncertainty": round(U.uncertainty(p.get("prob_sensitive"), m.get("auroc"),
                                                  p.get("n_train")), 4)}
        if neighbours and mod.nn is not None and drug in mod.nn["auc"].columns:
            v = mod.nn["auc"].loc[neighbours["specimens"], drug].dropna()
            if len(v) >= 5:
                allv = mod.nn["auc"][drug].dropna()
                row["neighbour_measured_auc_median"] = round(float(v.median()), 2)
                row["neighbour_n"] = int(len(v))
                row["cohort_measured_auc_median"] = round(float(allv.median()), 2)
                row["neighbour_percentile_in_cohort"] = round(float((allv <= v.median()).mean()), 3)
        rows.append(row)
    rows.sort(key=lambda r: -(r["prob_sensitive"] or 0))
    return _res("drug_response", "honest_cv", "rna_derived",
                {"n_drugs": len(rows), "top": rows[:top], "all": rows},
                {"note": "predicted ex-vivo response; not a predicted clinical benefit"})


# ------------------------------------------------ 2. cell-state coverage ----
def cell_state_agent(state_result, top=15):
    if not state_result or not state_result.get("per_drug"):
        return _res("cell_state_coverage", "descriptive_aggregate", "rna_derived",
                    {"available": False,
                     "reason": "no per-cell-state data for this sample (bulk input, or no state passed "
                               "the cell-count floor)"})
    rows = []
    for drug, v in state_result["per_drug"].items():
        rows.append({"inhibitor": drug, "coverage_blast": v["coverage_blast"],
                     "coverage_LSC_like": v["coverage_LSC_like"], "worst_state": v["worst_state"],
                     "worst_state_sens": v["worst_state_sens"], "dispersion": v["dispersion"],
                     "bulk_vs_sc": v["bulk_vs_sc"], "prob_sensitive_bulk": v["prob_sensitive_bulk"]})
    rows.sort(key=lambda r: -((r["coverage_blast"] or 0) + (r["coverage_LSC_like"] or 0)))
    disagree = sorted(rows, key=lambda r: -abs(r["bulk_vs_sc"] or 0))[:5]
    return _res("cell_state_coverage", "descriptive_aggregate", "rna_derived",
                {"available": True, "states": state_result["states"],
                 "compartment_note": state_result["compartment_note"],
                 "top": rows[:top], "largest_bulk_vs_singlecell_disagreement": disagree},
                {"note": "compartment coverage, not a malignant-cell call: no per-cell genotype is used"})


# --------------------------------------------------- 3. molecular mechanism ----
def mechanism_agent(mech_by_drug, top=15):
    rows = [{"inhibitor": d, "mechanistic_score": v.get("mechanistic_score"),
             "target_expression_pct": v.get("target_expression_pct"),
             "pathway_readout_pct": v.get("pathway_readout_pct"),
             "bcl2_dependency": v.get("bcl2_dependency"),
             "genetic_activation": v.get("genetic_activation"),
             "resistance_flags": v.get("resistance_flags")}
            for d, v in mech_by_drug.items()]
    rows.sort(key=lambda r: -(r["mechanistic_score"] or 0))
    return _res("molecular_mechanism", "classifier_call", "rna_derived",
                {"n_drugs": len(rows), "top": rows[:top], "all": rows},
                {"note": "independent of the empirical response model by design; agreement between the "
                         "two is the strong evidence, disagreement is informative"})


# --------------------------------------------------------- 4. pharmacology ----
def pharmacology_agent(curves, drugs, preds=None):
    """Does the measured sensitivity live at concentrations a patient could plausibly reach, and is it
    on-target? The curation answers the first half; the curve fits answer the second."""
    rows = []
    for d in drugs:
        ann = TG.get(d)
        rel = TG.assay_reliance(curves, d) or {}
        concern = []
        if ann["exposure"] == "not_established":
            concern.append("no established human exposure for this compound")
        elif ann["exposure"] == "borderline":
            concern.append("clinical exposure only borderline overlaps the assay window")
        f = rel.get("frac_sensitive_ic50_top_decade")
        if f is not None and f >= 0.5:
            concern.append("sensitivity mostly appears only in the top decade of the tested "
                           "concentration range (%.0f%% of sensitive specimens) -- selectivity doubtful" % (100 * f))
        if len(ann["targets"]) >= 6:
            concern.append("multi-kinase compound (%d annotated targets): a response cannot be "
                           "attributed to any one target" % len(ann["targets"]))
        rows.append({"inhibitor": d, "exposure": ann["exposure"],
                     "exposure_score": TG.exposure_score(d), "n_targets": len(ann["targets"]),
                     "assay_window_uM": [rel.get("min_conc_uM"), rel.get("max_conc_uM")],
                     "frac_sensitive_ic50_top_decade": f, "concerns": concern})
    rows.sort(key=lambda r: (len(r["concerns"]), -r["exposure_score"]))
    return _res("pharmacology", "descriptive_aggregate", "independent",
                {"n_drugs": len(rows), "clean": [r for r in rows if not r["concerns"]][:15],
                 "flagged": [r for r in rows if r["concerns"]][:15], "all": rows})


# ----------------------------------------------------- 5. clinical evidence ----
def clinical_evidence_agent(drugs, clinical=None):
    tiers = {t: [] for t in TG.TIERS}
    for d in drugs:
        a = TG.get(d)
        tiers[a["clinical_tier"]].append(
            {"inhibitor": d, "mechanism": a["mechanism"], "analogue": a["analogue"],
             "infeasibility": U.infeasibility(d, clinical)[0]})
    return _res("clinical_evidence", "deterministic_fact", "independent",
                {"tiers": {t: {"label": TG.TIER_LABEL[t], "n": len(v), "drugs": v} for t, v in tiers.items()},
                 "eligibility_assessed": bool(clinical)},
                {"note": "rankings are produced per tier; a research-only compound is never presented "
                         "as a treatment option"})


# --------------------------------------------------------- 6. combination ----
def combination_agent(state_result, scored, max_pairs=8, min_gain=0.15):
    """Complementary-coverage hypotheses ONLY.

    BeatAML measured single agents. Adding two single-agent sensitivity scores is not evidence of
    combination efficacy and this agent refuses to do it. What it can legitimately propose is a pair
    where the cell states one drug fails to cover are covered by the other, and the two act through
    different pathways -- a resistance-escape hypothesis for laboratory testing.
    """
    if not state_result or not state_result.get("per_state"):
        return _res("combination", "descriptive_aggregate", "rna_derived",
                    {"available": False, "reason": "combination reasoning requires per-cell-state "
                                                   "response, which needs single-cell input"})
    info = {s["state"]: s for s in state_result["states"]}
    blast = [s for s, v in info.items() if v["blast_compartment"] and v["fraction"] >= SM.MIN_FRAC]
    if len(blast) < 3:
        return _res("combination", "descriptive_aggregate", "rna_derived",
                    {"available": False, "reason": "too few blast-compartment states to reason about coverage"})
    w = np.array([info[s]["fraction"] for s in blast]); w = w / w.sum()

    cand = [d for d, v in scored.items()
            if TG.get(d)["clinical_tier"] in ("approved_AML", "approved_other")
            and v["utility"] > 0][:25]
    cov = {}
    for d in cand:
        ps = state_result["per_state"].get(d) or {}
        cov[d] = np.array([1.0 if (ps.get(s) is not None and ps[s] > 0) else 0.0 for s in blast])
    pairs = []
    for i, a in enumerate(cand):
        for b in cand[i + 1:]:
            if TG.get(a)["family_group"] == TG.get(b)["family_group"]:
                continue                       # same pathway: not a complementary-coverage hypothesis
            ca, cb = float(w @ cov[a]), float(w @ cov[b])
            cu = float(w @ np.maximum(cov[a], cov[b]))
            gain = cu - max(ca, cb)
            if gain >= min_gain:
                pairs.append({"pair": [a, b], "coverage_a": round(ca, 3), "coverage_b": round(cb, 3),
                              "coverage_union": round(cu, 3), "gain": round(gain, 3),
                              "pathways": [TG.get(a)["family_group"], TG.get(b)["family_group"]],
                              "states_rescued": [s for s, x, y in zip(blast, cov[a], cov[b])
                                                 if x == 0 and y == 1 or y == 0 and x == 1]})
    pairs.sort(key=lambda p: -p["gain"])
    return _res("combination", "descriptive_aggregate", "rna_derived",
                {"available": True, "n_pairs": len(pairs), "pairs": pairs[:max_pairs],
                 "basis": "complementary cell-state coverage across different target pathways"},
                {"refusal": "single-agent sensitivity scores were NOT added; BeatAML2 contains no "
                            "combination measurements, so no synergy is claimed",
                 "status": "hypotheses for laboratory testing, not treatment recommendations"})


# ------------------------------------------------------------- 7. skeptic ----
def skeptic_agent(mod, drug, pred, mech, state_metrics, state_result, ood_q, drug_tier,
                  patient_state_axis_q=None, curves=None):
    """Nine specific challenges. Each returns a concrete objection or nothing."""
    ch = []
    oof = (getattr(mod, "oof_metrics", {}) or {}).get(drug, {})
    gm = mod.group_models.get(mod.drug_group.get(drug), {})
    wstate = (gm.get("w") or {}).get("state", 0.0)

    if wstate >= 0.40 and patient_state_axis_q is not None and (patient_state_axis_q <= 0.1 or patient_state_axis_q >= 0.9):
        ch.append({"challenge": "differentiation_state_confound", "severity": "high",
                   "detail": "this drug's shared model puts %.0f%% of its weight on the differentiation-"
                             "state block and this patient sits at the %.0fth percentile of that axis; "
                             "the prediction may be a differentiation-state effect rather than a "
                             "drug-specific one" % (100 * wstate, 100 * patient_state_axis_q)})
    if (pred.get("n_train") or 0) < 200:
        ch.append({"challenge": "small_training_set", "severity": "medium",
                   "detail": "only %d BeatAML specimens back this inhibitor" % (pred.get("n_train") or 0)})
    if drug_tier == "wave_conditional":
        ch.append({"challenge": "assay_batch", "severity": "high",
                   "detail": "this inhibitor's AUC distribution shifted between BeatAML acquisition "
                             "waves; only leave-wave-out validation is meaningful for it"})
    ms = (mech or {}).get("mechanistic_score")
    ps = pred.get("prob_sensitive")
    if ms is not None and ps is not None and abs(ms - ps) >= 0.4:
        ch.append({"challenge": "modality_conflict", "severity": "medium",
                   "detail": "empirical P(sensitive)=%.2f but mechanistic score=%.2f; the two lines of "
                             "evidence disagree" % (ps, ms)})
    if curves is not None:
        rel = TG.assay_reliance(curves, drug) or {}
        f = rel.get("frac_sensitive_ic50_top_decade")
        if f is not None and f >= 0.5:
            ch.append({"challenge": "curve_quality", "severity": "medium",
                       "detail": "%.0f%% of sensitive specimens only reach IC50 in the top decade of "
                                 "the tested range" % (100 * f)})
    if ood_q is not None and ood_q >= 0.9:
        ch.append({"challenge": "out_of_distribution", "severity": "high",
                   "detail": "this patient sits at the %.0fth percentile of distance from the BeatAML "
                             "training distribution; the prediction is an extrapolation" % (100 * ood_q)})
    ann = TG.get(drug)
    if ann["clinical_tier"] == "research" and not ann["analogue"]:
        ch.append({"challenge": "no_clinical_route", "severity": "high",
                   "detail": "research-only compound with no clinically available analogue"})
    if (oof.get("auroc") or 0) < 0.65:
        ch.append({"challenge": "weak_model", "severity": "high",
                   "detail": "the model's own held-out AUROC for this inhibitor is %.2f" % (oof.get("auroc") or 0)})
    # on-target effect in normal haematopoiesis
    if state_result and state_result.get("per_state", {}).get(drug):
        ps_map = state_result["per_state"][drug]
        info = {s["state"]: s for s in state_result["states"]}
        norm = [v for s, v in ps_map.items() if s in info and not info[s]["blast_compartment"]]
        blast = [v for s, v in ps_map.items() if s in info and info[s]["blast_compartment"]]
        if len(norm) >= 3 and len(blast) >= 3 and np.mean(norm) >= np.mean(blast):
            ch.append({"challenge": "normal_haematopoiesis", "severity": "medium",
                       "detail": "predicted sensitivity is no lower in the presumed-normal compartment "
                                 "(%.2f) than in the blast compartment (%.2f): on-target myelo-"
                                 "suppression is plausible" % (np.mean(norm), np.mean(blast))})
    return ch


# ----------------------------------------------------------- 8. reporting ----
def reporting_agent(ranked, per_drug_evidence, abstained, patient_meta=None, top_per_tier=5):
    """A structured recommendation. No free narrative: every field is traceable to an evidence item."""
    out = {}
    for tier, blk in ranked.items():
        items = []
        for r in blk["ranked"][:top_per_tier]:
            d = r["inhibitor"]
            ev = per_drug_evidence.get(d, {})
            items.append({
                "inhibitor": d, "utility": r["utility"], "components": r["components"],
                "predicted_percentile": (ev.get("response") or {}).get("percentile"),
                "calibrated_probability": (ev.get("response") or {}).get("prob_sensitive"),
                "model_reliability_auroc": (ev.get("response") or {}).get("model_oof_auroc"),
                "beataml_support": {"n_train": (ev.get("response") or {}).get("n_train"),
                                    "neighbour_measured_auc_median":
                                        (ev.get("response") or {}).get("neighbour_measured_auc_median"),
                                    "cohort_measured_auc_median":
                                        (ev.get("response") or {}).get("cohort_measured_auc_median")},
                "coverage": {"blast": (ev.get("state") or {}).get("coverage_blast"),
                             "LSC_like": (ev.get("state") or {}).get("coverage_LSC_like"),
                             "escape_state": (ev.get("state") or {}).get("worst_state")},
                "mechanism": {"score": (ev.get("mechanism") or {}).get("mechanistic_score"),
                              "target_expression_pct": (ev.get("mechanism") or {}).get("target_expression_pct"),
                              "genetic_activation": (ev.get("mechanism") or {}).get("genetic_activation")},
                "clinical_status": TG.TIER_LABEL[TG.get(d)["clinical_tier"]],
                "challenges": ev.get("challenges") or [],
                "interpretation": ("candidate for trial matching or laboratory validation, "
                                   "not direct treatment selection"),
            })
        out[tier] = {"label": blk["label"], "n_considered": blk["n"], "recommendations": items}
    return _res("reporting", "descriptive_aggregate", "independent",
                {"by_tier": out, "abstained": abstained,
                 "patient": patient_meta or {},
                 "global_caveat": ("ex-vivo sensitivity is an experimentally grounded prioritisation "
                                   "signal, not an estimate of clinical benefit: culture conditions, "
                                   "pharmacokinetics, the marrow microenvironment, combination therapy "
                                   "and toxicity are not represented")})
