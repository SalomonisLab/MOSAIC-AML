#!/usr/bin/env python3
"""The survival hook the patient pipeline calls: expression (+ whatever clinical data exists) in,
a survival curve out.

Two things this module refuses to do, both of which would make the output look better than it is:

1. **It will not silently use the strong model without the data the strong model needs.** `full`
   (C-index 0.752 on the sealed hold-out) requires age and ELN 2017. Most uploads have neither. When
   they are missing the layer falls back to `molecular` (0.716) and *says so*, rather than imputing a
   median age and quoting the strong model's accuracy.
2. **It will not report a single number of months as the answer.** Group-level median survival is
   accurate to about six weeks, but actual survival inside one predicted-risk band spans roughly a
   year. The per-patient median is returned, labelled, next to the spread that makes it honest.

Kept behind a try/except at the call site like the drug layer: a missing model must cost the survival
section, never the mutation panel.
"""
from __future__ import annotations
import os, sys, json, pickle
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
MODEL = os.path.join(HERE, "survival_model.pkl")


def available():
    return os.path.exists(MODEL)


def _load():
    import json
    with open(MODEL, "rb") as f:
        b = pickle.load(f)
    tv = os.path.join(os.path.dirname(HERE), "deliverables", "survival_time_validation.json")
    b["card"]["_time"] = json.load(open(tv)) if os.path.exists(tv) else {}
    return b



def align_and_z(bundle, expr, sym2ens=None, cohort="sc"):
    """Map a gene-keyed expression vector onto the model's ENSG space and z-score it.

    Returns (Z, n_matched). Owning the alignment here is deliberate: the first smoke test silently fed
    all-zero vectors because the atlas is keyed by gene SYMBOL and the model by ENSG, and an all-zero
    input still came back as a confident 97th-percentile high-risk prognosis. Anything that can quietly
    turn a gene-ID mismatch into a confident prediction about how long someone will live has to be
    caught here, not trusted to the caller.
    """
    import pandas as pd
    fs = bundle["feature_space"]
    gpos = {g: i for i, g in enumerate(fs.genes)}
    s2e = sym2ens or bundle.get("sym2ens") or {}
    vals = expr.values if hasattr(expr, "values") else np.asarray(expr)
    keys = list(expr.index if hasattr(expr, "index") else expr.columns)
    A = np.zeros((1, len(fs.genes)))
    n = 0
    for k, v in zip(keys, np.asarray(vals).ravel()):
        k = str(k)
        j = gpos.get(k if k in gpos else s2e.get(k))
        if j is not None:
            A[0, j] = float(v); n += 1
    ref = cohort if cohort in fs.ref else "beataml"
    return fs.z(A, ref), n


# Evidence-based, not a guess: withholding genes from one sample drifts its risk percentile 0.06 ->
# 0.08 (90% coverage) -> 0.23 (75%) -> 0.58 (50%), i.e. below ~80% the risk GROUP starts to flip. Refuse
# under 0.80 and warn between 0.80 and 0.95.
MIN_GENE_FRACTION = 0.80
WARN_GENE_FRACTION = 0.95


def _blocks_for(bundle, z_row, mutations=None, clinical=None):
    """Build the feature blocks for ONE patient from an already-z-scored expression row."""
    from amlmm.survival import data as SD
    import pandas as pd
    fs = bundle["feature_space"]
    P = fs.pca.transform(z_row[:, fs.sel])
    S, _ = fs._state_block(z_row)
    # place the mutation caller's output into the exact columns the model was fitted on; anything the
    # caller did not report stays 0 (= not observed present), and the trailing column records how much
    # of the panel was unreported so the model can discount a thin call set
    cols = bundle.get("mut_columns") or []
    M = None
    if cols:
        M = np.zeros((1, len(cols)))
        drivers = cols[:-1]
        if mutations:
            hit = 0
            for j, name in enumerate(drivers):
                v = mutations.get(name)
                if v is not None:
                    M[0, j] = float(v); hit += 1
            M[0, -1] = 1.0 - hit / max(1, len(drivers))
        else:
            M[0, -1] = 1.0
    out = {"rna": P, "state": S}
    if M is not None:
        out["mut"] = M
    if clinical:
        cl = pd.DataFrame([{"ageAtDiagnosis": clinical.get("age"),
                            "ELN2017": clinical.get("eln"),
                            "consensus_sex": clinical.get("sex"),
                            "%.Blasts.in.BM": clinical.get("blasts_bm"),
                            "%.Blasts.in.PB": clinical.get("blasts_pb"),
                            "wbcCount": clinical.get("wbc"),
                            "isRelapse": clinical.get("is_relapse"),
                            "isDenovo": clinical.get("is_denovo"),
                            "isTransformed": clinical.get("is_transformed"),
                            "priorMDS": clinical.get("prior_mds")}])
        out["clin"], _ = SD.clinical_block(cl)
        out["age_eln"], _ = SD.age_eln_block(cl)
        age = out["age_eln"][:, [0]]
        kn = bundle.get("age_knots") or [45.0, 60.0, 71.0]
        out["age_spline"] = np.hstack([age] + [np.clip(age - q, 0, None) ** 3 for q in kn])
        ty = str(clinical.get("induction") or "")
        out["txbase"] = np.array([[float("chemo" in ty.lower()),
                                   float("palliat" in ty.lower()),
                                   float(not ty)]])
    return out


_TCGA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "deliverables", "validation_tcga_laml.json")


def _external_validation(arm):
    """What the SAME arm scored on a cohort the model was never fitted on.

    A hold-out C-index says the model did not memorise these patients. It does not say the model did not
    memorise this *cohort*, and a report that quotes only the internal number invites exactly that
    over-reading. Carrying the TCGA-LAML transfer next to it makes the generalisation gap visible on the
    patient's own report instead of living in a methods file nobody opens. Absent file -> absent key,
    never a fabricated number.
    """
    try:
        with open(_TCGA, encoding="utf-8") as f:
            v = json.load(f)
        a = (v.get("arms") or {}).get(arm)
        if not a or a.get("c_index") is None:
            return None
        g = (v.get("gain_over_age_eln") or {}).get(arm) or {}
        return {"cohort": v.get("cohort"), "n": v.get("n"), "deaths": v.get("deaths"),
                "c_index": a.get("c_index"), "c_index_ci95": a.get("c_index_ci95"),
                "auc_2y": a.get("auc_2y"),
                "gain_over_age_cytogenetics": g.get("delta_c"),
                "note": ("coefficients frozen; only the per-gene z-reference was cohort-matched. The "
                         "drop from the internal hold-out is the cost of transferring across platforms "
                         "and treatment eras, and is the more realistic expectation for a new cohort.")}
    except Exception:
        return None


def _rmst(surv_t, surv_S, horizon):
    """Restricted mean survival time: the area under the survival curve up to `horizon`.

    Preferred over the median because it is DEFINED even when the curve never reaches 0.5 — the
    low-risk group's median was previously reported as 6.9 years, extrapolated past the end of
    follow-up. RMST says "on average, this many of the next N years survived", which is both
    well-defined and easier to act on.
    """
    import numpy as _np
    g = _np.concatenate([[0.0], surv_t[surv_t <= horizon], [horizon]])
    s = _np.concatenate([[1.0], surv_S[surv_t <= horizon],
                         [_np.interp(horizon, surv_t, surv_S, left=1.0,
                                     right=surv_S[-1] if len(surv_S) else 1.0)]])
    return float(_np.sum(_np.diff(g) * s[:-1]))


def _pick_arm(bundle, blocks, clinical):
    """The best arm whose inputs are actually present."""
    have_clin = bool(clinical) and clinical.get("age") is not None
    for arm in (("deployed", "full", "molecular") if have_clin else ("molecular",)):
        if arm in bundle["models"] and all(b in blocks for b in bundle["arm_blocks"][arm]):
            return arm
    return next(iter(bundle["models"]))


def _pick_stratum(bundle, clinical):
    """Which treatment stratum this patient belongs to, or None to use the pooled model.

    A model fitted on the pooled cohort and applied to non-intensively-treated patients scores C-index
    0.554; fitted WITHIN that stratum it scores 0.681. So when the caller states the induction type we
    use the matching stratum model, and when they do not we stay pooled -- guessing the stratum would be
    worse than not stratifying, and most uploads carry no treatment information at all.
    """
    strata = bundle.get("strata") or {}
    if not strata or not clinical:
        return None
    ty = str(clinical.get("induction") or "").strip().lower()
    if not ty:
        return None
    name = "intensive" if ("chemo" in ty or "intensive" in ty) and "non" not in ty else "non_intensive"
    if name not in strata:
        return None
    st = strata[name]
    return {"name": name, "model": st["model"], "risk_ref": st["risk_ref"],
            "n_train": st.get("n_train"), "events_train": st.get("events_train")}


def run_for_expression(z_row, mutations=None, clinical=None, cohort="beataml",
                       n_genes_matched=None, specimen_class=None):
    """z_row: (1 x genes) expression, z-scored against the reference named by `cohort`
    ("beataml" for a bulk specimen, "sc" for a single-cell bulk-equivalent)."""
    if not available():
        return {"available": False, "reason": "survival_model.pkl not built (run train_survival_model.py)"}
    b = _load()
    card = b["card"]
    ngen = len(b["feature_space"].genes)
    if n_genes_matched is not None and n_genes_matched < MIN_GENE_FRACTION * ngen:
        return {"available": False,
                "reason": ("only %d of the model's %d genes matched this sample (%.0f%%, floor %.0f%%) — "
                           "most likely a gene-identifier mismatch. Refusing to report a prognosis from "
                           "a mostly-empty expression vector."
                           % (n_genes_matched, ngen, 100 * n_genes_matched / ngen,
                              100 * MIN_GENE_FRACTION))}
    if not np.isfinite(z_row).all() or float(np.nanstd(z_row)) < 1e-9:
        return {"available": False, "reason": "the expression vector is empty or constant"}
    # A survival model fitted on AML patients has nothing to say about a healthy donor. Without this
    # guard a healthy pooled-CD34 control scored at the 97th risk percentile with a 0.4% one-year
    # survival, which is not a cautious answer — it is a wrong one.
    if str(specimen_class).lower() == "control":
        return {"available": False,
                "reason": ("the control gate called this specimen healthy; the survival model is fitted "
                           "on AML patients and does not apply to a non-leukaemic sample")}
    blocks = _blocks_for(b, z_row, mutations, clinical)
    arm = _pick_arm(b, blocks, clinical)
    m = b["models"][arm]
    names = b["arm_blocks"][arm]
    sub = {k: blocks[k] for k in names if k in blocks}
    if len(sub) != len(names):
        return {"available": False, "reason": "no arm could be fed with the data supplied"}

    # Treatment stratum, when the caller actually stated the induction type. Only the `deployed` arm has
    # stratum-specific models, and only they are used -- a stratum model applied through a different
    # arm's block list would be a silent mismatch.
    strat = _pick_stratum(b, clinical) if arm == "deployed" else None
    if strat is not None:
        m = strat["model"]

    raw = float((m.risk(sub) if hasattr(m, "stack_blocks") else m.risk(sub[names[0]]))[0])
    ref = strat["risk_ref"] if strat is not None else b["risk_ref"].get(arm)
    # Cohort-matched risk: for single-cell input, rank among single-cell samples and read off the
    # BeatAML risk at that rank. The baseline hazard is BeatAML's — it is the only one tied to observed
    # survival — so the score fed into it has to be on BeatAML's scale, not on an sc sample's.
    ref_sc = (b.get("risk_ref_sc") or {}).get(arm) if cohort == "sc" else None
    if ref_sc is not None and len(ref_sc) and ref is not None and len(ref):
        pct = float(np.searchsorted(ref_sc, raw, side="right") / len(ref_sc))
        risk = float(np.quantile(ref, min(max(pct, 1e-3), 1 - 1e-3)))
    else:
        risk = raw
        pct = None if ref is None else float(np.searchsorted(ref, risk, side="right") / len(ref))
    H = b["horizons"]
    S = (m.survival(sub, H) if hasattr(m, "stack_blocks") else m._final.survival(
        np.array([[risk]]), H))[0]
    med = (m.median_survival(sub) if hasattr(m, "stack_blocks")
           else m._final.median_survival(np.array([[risk]])))[0]

    # RMST at each horizon, from this patient's own predicted curve
    grid = np.linspace(0.05, max(H), 200)
    Sg = (m.survival(sub, grid) if hasattr(m, "stack_blocks")
          else m._final.survival(np.array([[risk]]), grid))[0]
    rmst_out = {"%gy" % h: round(_rmst(grid, Sg, h), 2) for h in H}

    tert = "intermediate"
    if ref is not None and len(ref):
        lo, hi = np.quantile(ref, [1 / 3, 2 / 3])
        tert = "low" if risk < lo else ("high" if risk >= hi else "intermediate")
    tv = card.get("_time") or {}
    spread = (tv.get("individual_spread_within_risk_group") or {}).get(tert)

    ho = card["holdout"].get(arm, {})
    return {
        "available": True, "arm": arm,
        "arm_note": ("age and ELN 2017 were supplied, so the stronger combined model was used"
                     if arm == "full" else
                     "no age / ELN 2017 supplied, so the molecular-only model was used — it is weaker "
                     "(hold-out C-index %.3f vs %.3f for the combined model)"
                     % (card["holdout"].get("molecular", {}).get("c_index", float("nan")),
                        card["holdout"].get("full", {}).get("c_index", float("nan")))),
        "risk_score": round(risk, 4), "risk_score_raw": round(raw, 4), "cohort_reference": cohort,
        "genes_matched": n_genes_matched,
        "gene_coverage_warning": (None if (n_genes_matched is None
                                           or n_genes_matched >= WARN_GENE_FRACTION * ngen)
                                  else "only %.0f%% of the model's genes matched; the risk score drifts "
                                       "as coverage falls, so treat the risk group as approximate"
                                       % (100 * n_genes_matched / ngen)),
        "risk_percentile": None if pct is None else round(pct, 3),
        "risk_group": tert,
        "survival_probability": {"%gy" % h: round(float(s), 3) for h, s in zip(H, S)},
        "rmst_years": rmst_out,
        "rmst_note": ("restricted mean survival — the average number of the next N years survived by "
                      "comparable patients. Defined even when the median is never reached, which is "
                      "why it is reported ahead of the median."),
        "median_survival_years": None if not np.isfinite(med) else round(float(med), 2),
        "median_survival_note": ("not reached inside the training follow-up — more than half of "
                                 "comparable patients were still alive" if not np.isfinite(med) else
                                 "a point estimate; see the spread below before quoting it"),
        "observed_spread_in_this_risk_group": spread,
        "model_accuracy": {"c_index_holdout": ho.get("c_index"),
                           "auc_2y_holdout": ho.get("auc_2y"),
                           "c_index_gain_over_age_eln": (card["incremental"].get(arm) or {}).get("delta_c"),
                           "n_train": card["cohort"]["final_patients"],
                           "events_train": card["cohort"]["events"],
                           "external_validation": _external_validation(arm),
                           "treatment_stratum": (None if strat is None else {
                               "stratum": strat["name"], "n_train": strat["n_train"],
                               "events_train": strat["events_train"],
                               "cv_c_index": (0.729 if strat["name"] == "intensive" else 0.681),
                               "why": ("fitted within this treatment stratum. The pooled model applied "
                                       "to non-intensively-treated patients scores 0.554; fitted within "
                                       "that stratum it scores 0.681. ELN 2022 scores 0.462 there, at "
                                       "or below chance.")}),
                           "stratum_note": (None if strat is not None else
                                            "no induction type supplied, so the pooled model was used; "
                                            "supplying it selects a treatment-stratified model")},
        "caveat": ("Trained on BeatAML2 (bulk RNA, initial-diagnosis specimens). This is a prognostic "
                   "estimate from data available at diagnosis, not a statement about what will happen "
                   "to this patient, and it does not account for the treatment they go on to receive."),
        # A gate result only silences this warning if it actually SAYS something. `not_excluded` does
        # not: the gate calls ~50% of genuinely healthy specimens positive, so treating its positive
        # call as confirmation that this is a patient would launder a coin flip into a premise.
        "assumes_aml": (None if str(specimen_class).lower() == "control" else
                        (("the control gate did not exclude disease, but at its operating point it "
                          "calls about half of healthy specimens positive, so it has NOT established "
                          "that this is an AML patient — this number still assumes it.")
                         if specimen_class else
                         ("no healthy-vs-diseased gate ran for this input (that needs single cells), so "
                          "this number ASSUMES the sample is from an AML patient. Scored against a "
                          "healthy donor it will still return a risk percentile, and that percentile "
                          "will be meaningless."))),
    }



# --------------------------------------------------- pipeline entry points ----
def _mut_probs(mut_preds, observed):
    """The mutation caller's output as {category: probability}, with observed genotypes pinned to 1."""
    out = {}
    for x in (mut_preds or []):
        c = x.get("category") or x.get("mutation")
        if c:
            out[str(c)] = float(x.get("probability") or 0.0)
    for m in (observed or []):
        out[str(m)] = 1.0
    return out or None


def run_for_sample(expr_series, cohort="sc", mut_preds=None, observed=None, clinical=None,
                   specimen_class=None):
    """One gene-keyed expression vector (single-cell bulk-equivalent or bulk RNA) -> survival block."""
    if not available():
        return {"available": False, "reason": "survival_model.pkl not built"}
    b = _load()
    Z, n = align_and_z(b, expr_series, cohort=cohort)
    return run_for_expression(Z, mutations=_mut_probs(mut_preds, observed), clinical=clinical,
                              cohort=cohort, n_genes_matched=n, specimen_class=specimen_class)


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser(description="Run the survival layer on one atlas sample (smoke test).")
    ap.add_argument("--atlas-sample", default="CCHMC::1009_AfInv16_29M")
    ap.add_argument("--age", type=float)
    ap.add_argument("--eln", choices=["Favorable", "Intermediate", "Adverse"])
    a = ap.parse_args()
    import pandas as pd
    import predict_drugs as PD
    from amlmm.drug import statemodel as SM

    counts, _ = PD.atlas_sample(a.atlas_sample)
    lin = SM.cp10k(counts.values.sum(0, keepdims=True))[0]
    ser = pd.Series(lin, index=[str(c) for c in counts.columns])
    b = _load()
    Z, n = align_and_z(b, ser, cohort="sc")
    clin = {"age": a.age, "eln": a.eln} if a.age is not None else None
    print(json.dumps(run_for_expression(Z, clinical=clin, cohort="sc", n_genes_matched=n), indent=1))
