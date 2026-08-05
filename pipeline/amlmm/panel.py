"""Per-modality expert agent panel + cross-modality arbiter (MOSAIC-AML core).

A "tumor board": each witness owns one kind of evidence and returns a structured,
grounded read; a chair (arbiter) reconciles them into a consensus that accounts for
every witness. Three witness kinds:

  * predictive  -- a feature block (composition / a modality). Grounded in that
                   block's HONEST donor-grouped CV (balanced acc + permutation p +
                   top features). Measured blocks weighted above imputed-from-RNA.
  * genetic     -- mutation/cytogenetic context + targetable drivers (NOT a subtype
                   predictor; provides risk + therapy priors). An independent axis.
  * udon        -- conserved UDON programs (control-normalized) and their subtype
                   associations; flags batch-dominated programs. Discovery evidence.

Every LLM call validates + falls back to a deterministic read, so a gateway outage
degrades rather than breaks. The arbiter is told to treat imputed-from-RNA modalities
as corroborating, NOT independent (shared RNA source).
"""
from __future__ import annotations
import json

import numpy as np

from .step import run_step, get
from .llm import LLMClient, LLMError
from . import steps as _steps   # noqa: F401  -- registers feasibility/assemble/classify
from . import genetics, udon, dataio

MEASURED = {"composition", "RNA"}
OPINION_KEYS = ("confidence", "reliability_weight", "summary", "caveats")
CONSENSUS_KEYS = ("consensus", "overall_confidence", "per_witness_consistency",
                  "targetable_summary", "recommended_validations", "conflicts", "caveat")


# ---------------------------------------------------------------- predictive witness
def _kind(block):
    return "measured" if block in MEASURED else "imputed_from_RNA"


def top_features(X, y, k=8):
    try:
        from sklearn.ensemble import RandomForestClassifier
        rf = RandomForestClassifier(n_estimators=200, class_weight="balanced_subsample",
                                    random_state=0, n_jobs=1)
        rf.fit(np.asarray(X.values, float), np.asarray(y).astype(str))
        order = np.argsort(rf.feature_importances_)[::-1][:k]
        cols = list(X.columns)
        return [str(cols[i]).split("::", 1)[-1] for i in order]
    except Exception:
        return []


def evidence_predictive(ctx, target, block, strategy="donor_kfold", permutations=40):
    a = run_step(get("assemble_features"), ctx,
                 {"target": target, "blocks": [block], "min_class_n": 8})
    if a.status != "ok":
        return {"witness": block, "kind": _kind(block), "status": a.status,
                "reason": a.metrics.get("reason", "")}
    run_step(get("classify"), ctx, {"strategy": strategy, "outer_splits": 5,
                                    "inner_splits": 3, "n_permutations": permutations})
    cv = ctx.getart("cv_result") or {}
    X, y = ctx.getart("X"), ctx.getart("y")
    pval = cv.get("permutation_pvalue")
    return {"witness": block, "kind": _kind(block), "status": "ok",
            "balanced_accuracy": cv.get("balanced_accuracy"), "permutation_pvalue": pval,
            "above_chance": (pval is not None and pval < 0.05),
            "classes": cv.get("classes"), "per_class_f1": cv.get("per_class_f1"),
            "top_features": top_features(X, y) if X is not None else []}


def assess_predictive(client, target, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0,
                "summary": f"no usable evidence ({ev.get('status')}: {ev.get('reason','')})",
                "caveats": "witness skipped"}
    indep = ("directly measured" if ev["kind"] == "measured"
             else "DERIVED from RNA -- not independent evidence")
    prompt = (
        f"You are the {ev['witness']} expert on a tumor-board panel predicting AML {target} "
        f"({indep}). Honest donor-grouped CV evidence for predicting {target} from THIS block alone:\n"
        f"  balanced_accuracy={ev['balanced_accuracy']}, permutation p={ev['permutation_pvalue']} "
        f"(above chance={ev['above_chance']}), classes={ev['classes']}\n"
        f"  per-class F1={ev['per_class_f1']}\n  top features={ev['top_features']}\n"
        "Summarize how much this block contributes, which classes it best distinguishes, what the "
        "features suggest, your confidence, and caveats. If p>=0.05 it adds little alone. "
        "reliability_weight must reflect measured-vs-imputed AND above-chance. "
        'Return ONLY JSON {"confidence":0..1,"reliability_weight":0..1,"summary":str,"caveats":str}.')
    try:
        return client.chat_json(prompt, required=OPINION_KEYS)
    except (LLMError, Exception):
        w = (0.9 if ev["kind"] == "measured" else 0.5) * (1.0 if ev.get("above_chance") else 0.3)
        return {"confidence": 0.6 if ev.get("above_chance") else 0.2, "reliability_weight": round(w, 2),
                "summary": f"balanced_acc {ev.get('balanced_accuracy')}, above_chance={ev.get('above_chance')}",
                "caveats": "LLM unavailable; deterministic fallback"}


# ---------------------------------------------------------------- genetic witness
def evidence_genetic(ctx, target):
    g = genetics.summarize_genetics(ctx)
    g.update({"witness": "genetic", "kind": "independent", "status": "ok"})
    return g


def assess_genetic(client, target, ev):
    prompt = (
        "You are the genetic witness on an AML tumor-board panel. You provide RISK CONTEXT and "
        "TARGETABLE DRIVERS (therapy priors), NOT a subtype prediction. Cohort genetics "
        f"(from sparse karyotype + driver calls; {ev.get('n_with_genetic_data')}/{ev.get('n_samples')} "
        "have any genetic data):\n"
        f"  mutation/cytogenetic prevalence={ev.get('mutation_prevalence')}\n"
        f"  targetable drivers present (with therapy prior)={ev.get('targetable_present')}\n"
        f"  ELN risk distribution={ev.get('eln_distribution')}\n"
        "Summarize the genetic risk landscape and which targetable drivers are present (with their "
        "therapy priors). Genetics is an INDEPENDENT, actionable axis, but down-weight reliability by "
        "the low coverage. "
        'Return ONLY JSON {"confidence":0..1,"reliability_weight":0..1,"summary":str,"caveats":str}.')
    try:
        return client.chat_json(prompt, required=OPINION_KEYS)
    except (LLMError, Exception):
        cov = (ev.get("n_with_genetic_data") or 0) / max(ev.get("n_samples") or 1, 1)
        return {"confidence": round(0.5 * cov + 0.3, 2), "reliability_weight": round(0.6 * cov + 0.2, 2),
                "summary": f"targetable: {list(ev.get('targetable_present', {}))}; "
                           f"top mutations: {list(ev.get('mutation_prevalence', {}))[:5]}",
                "caveats": "LLM unavailable; deterministic fallback; sparse genetic coverage"}


# ---------------------------------------------------------------- UDON witness
def evidence_udon(ctx, target):
    u = udon.udon_associations(ctx)
    u.update({"witness": "cell-state/UDON", "kind": "discovery",
              "status": "ok" if u.get("available") else "skipped"})
    return u


def assess_udon(client, target, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0,
                "summary": f"no UDON programs ({ev.get('reason','')})", "caveats": "witness skipped"}
    prompt = (
        "You are the cell-state / UDON witness on an AML tumor-board panel. UDON programs are conserved, "
        f"control-normalized disease programs. There are {ev.get('n_programs')} programs. "
        f"program -> top associated subtype (n, dataset_dominated)={ev.get('program_subtype_top')}. "
        f"{ev.get('n_batch_dominated')} programs are dataset-dominated (batch -- DISTRUST those). "
        "Summarize which programs mark which subtypes (phenotype-genotype links), which to trust vs "
        "distrust as batch, and your confidence. "
        'Return ONLY JSON {"confidence":0..1,"reliability_weight":0..1,"summary":str,"caveats":str}.')
    try:
        return client.chat_json(prompt, required=OPINION_KEYS)
    except (LLMError, Exception):
        return {"confidence": 0.5, "reliability_weight": 0.5,
                "summary": f"{ev.get('n_programs')} programs; {ev.get('n_batch_dominated')} batch-dominated",
                "caveats": "LLM unavailable; deterministic fallback"}


# ---------------------------------------------------------------- arbiter
def reconcile(client, target, items):
    brief = []
    for it in items:
        e, o = it["evidence"], it["opinion"]
        b = {"witness": e.get("witness"), "kind": e.get("kind"),
             "confidence": o.get("confidence"), "reliability_weight": o.get("reliability_weight"),
             "summary": o.get("summary")}
        if e.get("kind") in ("measured", "imputed_from_RNA"):
            b.update({"balanced_accuracy": e.get("balanced_accuracy"),
                      "permutation_pvalue": e.get("permutation_pvalue")})
        if e.get("witness") == "genetic":
            b["targetable_present"] = e.get("targetable_present")
        if e.get("witness") == "cell-state/UDON":
            b["n_batch_dominated"] = e.get("n_batch_dominated")
        brief.append(b)
    prompt = (
        f"You are the chair reconciling an AML tumor-board panel for {target}. Each witness assessed a "
        "different evidence source.\n"
        f"Panel:\n{json.dumps(brief, default=str)[:3800]}\n"
        "Produce a consensus that ACCOUNTS FOR EACH witness. Rules: weight INDEPENDENT axes -- the "
        "genetic witness and MEASURED blocks (composition, RNA) -- and any block with permutation p<0.05; "
        "treat imputed-from-RNA blocks as CORROBORATING, not independent (several agreeing is NOT "
        "independent confirmation -- shared RNA source); DISTRUST batch-dominated UDON programs. State the "
        "best-supported biology, where witnesses agree vs conflict, a targetable_summary (drivers/programs "
        "worth acting on) and concrete recommended_validations (e.g. flow markers from ADT, targeted "
        "assays, ex vivo drug tests). Be candid about overall confidence given small n and shared sourcing. "
        'Return ONLY JSON {"consensus":str,"overall_confidence":"low"|"medium"|"high",'
        '"per_witness_consistency":{"<witness>":"agree"|"neutral"|"conflict"},"targetable_summary":str,'
        '"recommended_validations":str,"conflicts":str,"caveat":str}.')
    try:
        return client.chat_json(prompt, required=CONSENSUS_KEYS)
    except (LLMError, Exception) as e:
        return {"consensus": "LLM unavailable; deterministic fallback",
                "overall_confidence": "low",
                "per_witness_consistency": {it["evidence"].get("witness"): "neutral" for it in items},
                "targetable_summary": "", "recommended_validations": "",
                "conflicts": "", "caveat": str(e)[:200]}


def run_panel(ctx, target, blocks=None, with_genetic=True, with_udon=True,
              strategy="donor_kfold", permutations=40, client=None):
    client = client or LLMClient()
    blocks = blocks or ["composition"]
    items = []
    for b in blocks:
        ev = evidence_predictive(ctx, target, b, strategy=strategy, permutations=permutations)
        items.append({"evidence": ev, "opinion": assess_predictive(client, target, ev)})
    if with_genetic:
        ev = evidence_genetic(ctx, target)
        items.append({"evidence": ev, "opinion": assess_genetic(client, target, ev)})
    if with_udon:
        ev = evidence_udon(ctx, target)
        items.append({"evidence": ev, "opinion": assess_udon(client, target, ev)})
    consensus = reconcile(client, target, items)
    report = {"mode": "modality_panel", "target": target, "strategy": strategy,
              "provenance": ctx.getart("provenance"), "panel": items, "consensus": consensus}
    ctx.save_json(report, "panel_report.json")
    _write_md(ctx, report)
    return report


# ================================================================ per-patient mode
PATIENT_CONSENSUS_KEYS = ("subtype_call", "targetable_therapies", "recommended_validations",
                          "overall_confidence", "conflicts", "rationale")


def _patient_predictive(ctx, sample_key, block, strategy, permutations):
    ev = evidence_predictive(ctx, "subtype", block, strategy=strategy, permutations=permutations)
    if ev.get("status") != "ok":
        return {"witness": block, "kind": _kind(block), "status": ev.get("status"),
                "reason": ev.get("reason", "")}
    oof = (ctx.getart("cv_result") or {}).get("oof", {}).get(sample_key)
    out = {"witness": block, "kind": _kind(block), "status": "ok",
           "held_out": oof is not None,
           "patient_prediction": oof.get("pred") if oof else None,
           "patient_probability": oof.get("prob") if oof else None,
           "true_label": oof.get("true") if oof else None,
           "cohort_balanced_accuracy": ev.get("balanced_accuracy"),
           "cohort_permutation_p": ev.get("permutation_pvalue")}
    if block == "composition":
        comp = ctx.tables["composition"]
        if sample_key in comp.index:
            row = comp.loc[sample_key]
            tot = float(row.sum()) or 1.0
            # explicit (-fraction, name) tiebreak so the top-6 selection is environment-stable
            # (plain sort_values uses unstable quicksort -> tie-boundary can differ across builds).
            frac = sorted(((str(k), float(v) / tot) for k, v in row.items()),
                          key=lambda kv: (-kv[1], kv[0]))[:6]
            out["patient_top_cellstates"] = {k: round(v, 3) for k, v in frac}
    return out


def _patient_genetic(ctx, sample_key):
    # None-guard, NOT `or`: a cached non-empty DataFrame raises "truth value is
    # ambiguous", which run_agent would swallow -> empty genetic evidence -> the
    # deterministic anchor silently disappears for every patient after the first.
    M = ctx.tables.get("mutations")
    if M is None:
        M = genetics.build_mutation_matrix(ctx)
    base = {"witness": "genetic", "kind": "independent", "status": "ok"}
    if sample_key not in M.index:
        return {**base, "present": [], "targetable": {}, "eln": None}
    row = M.loc[sample_key]
    present = [c.replace("mut_", "").replace("cyto_", "")
               for c in M.filter(regex="^(mut_|cyto_)").columns if row.get(c) == 1.0]
    targetable = {g.replace("mut_", ""): t for g, t in genetics.TARGETABLE.items()
                  if g in M.columns and row.get(g) == 1.0}
    return {**base, "present": present, "targetable": targetable,
            "eln": (None if "ELN_risk" not in M.columns or pd_isna(row.get("ELN_risk"))
                    else str(row.get("ELN_risk")))}


def _patient_udon(ctx, sample_key):
    prog = ctx.tables.get("udon_programs")
    base = {"witness": "cell-state/UDON", "kind": "discovery"}
    if prog is None or not {"Dataset", "Sample", "final_program"}.issubset(prog.columns):
        return {**base, "status": "skipped"}
    assoc = udon.udon_associations(ctx).get("program_subtype_top", {})
    p = prog.copy()
    p["sample_key"] = [f"{d}::{x}" for d, x in zip(p["Dataset"], p["Sample"])]
    mine = p[p["sample_key"] == sample_key]
    counts = mine["final_program"].value_counts()
    # stable order (count desc, then name) so the evidence hash is reproducible
    # even when two programs tie on pseudobulk count.
    ordered = sorted(counts.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))[:8]
    active = [{"program": str(k), "n_pseudobulks": int(v),
               "marks": assoc.get(str(k), {}).get("top_subtype")} for k, v in ordered]
    return {**base, "status": "ok", "n_pseudobulks": int(len(mine)), "active_programs": active}


def pd_isna(v):
    try:
        import pandas as pd
        return pd.isna(v)
    except Exception:
        return v is None


def assess_patient_predictive(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0,
                "summary": f"no evidence ({ev.get('status')})", "caveats": "skipped"}
    cells = (f"This patient's dominant cell-states: {ev.get('patient_top_cellstates')}. "
             if ev.get("patient_top_cellstates") else "")
    prompt = (f"Per-PATIENT read from the {ev['witness']} witness ({_kind(ev['witness'])}). For THIS patient "
              f"(held-out={ev['held_out']}), the model predicts subtype {ev['patient_prediction']} at "
              f"probability {ev['patient_probability']}. Cohort reliability of this block for subtype: "
              f"balanced_accuracy={ev['cohort_balanced_accuracy']}, permutation p={ev['cohort_permutation_p']}. "
              f"{cells}Summarize what this witness concludes FOR THIS PATIENT and how much to trust it "
              "(measured vs imputed, held-out, cohort reliability). "
              'Return ONLY JSON {"confidence":0..1,"reliability_weight":0..1,"summary":str,"caveats":str}.')
    try:
        return client.chat_json(prompt, required=OPINION_KEYS)
    except (LLMError, Exception):
        w = 0.85 if _kind(ev["witness"]) == "measured" else 0.5
        return {"confidence": 0.5, "reliability_weight": w,
                "summary": f"predicts {ev.get('patient_prediction')} (prob {ev.get('patient_probability')})",
                "caveats": "LLM fallback"}


def assess_patient_genetic(client, ev):
    prompt = ("Per-PATIENT read from the genetic witness (independent axis). FOR THIS PATIENT: "
              f"mutations/cytogenetics present={ev.get('present')}; targetable drivers (with therapy)"
              f"={ev.get('targetable')}; ELN risk={ev.get('eln')}. Summarize the genetic risk and the "
              "actionable targets for THIS patient (or note if no genetic data is available). "
              'Return ONLY JSON {"confidence":0..1,"reliability_weight":0..1,"summary":str,"caveats":str}.')
    try:
        return client.chat_json(prompt, required=OPINION_KEYS)
    except (LLMError, Exception):
        return {"confidence": 0.6 if ev.get("present") else 0.2,
                "reliability_weight": 0.85 if ev.get("present") else 0.1,
                "summary": f"present: {ev.get('present')}; targetable: {list(ev.get('targetable', {}))}",
                "caveats": "LLM fallback"}


def assess_patient_udon(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0, "summary": "no UDON data", "caveats": "skipped"}
    prompt = ("Per-PATIENT read from the cell-state/UDON witness. FOR THIS PATIENT the active conserved "
              f"disease programs (with what each marks across the cohort) are {ev.get('active_programs')} "
              f"over {ev.get('n_pseudobulks')} pseudobulks. Summarize which programs this patient activates "
              "and what they imply biologically. "
              'Return ONLY JSON {"confidence":0..1,"reliability_weight":0..1,"summary":str,"caveats":str}.')
    try:
        return client.chat_json(prompt, required=OPINION_KEYS)
    except (LLMError, Exception):
        return {"confidence": 0.6, "reliability_weight": 0.7,
                "summary": f"programs: {[p['program'] for p in ev.get('active_programs', [])]}",
                "caveats": "LLM fallback"}


# ============================================================ Phase B descriptive witnesses
# These six witnesses are CORROBORATING / DISCOVERY context, never subtype voters: each uses
# a domain string the arbiter's vote/anchor branches do not read, and exposes (optionally)
# three additive keys the arbiter harvests generically: `therapy_biomarkers` (KB biomarker
# keys, OBSERVED-only), `validation_claims` (KB claim_type keys), `descriptive_context`
# (short narrative strings). All emitted lists are sorted and floats rounded so the ledger's
# deterministic_evidence_hash is stable across runs.
ACTIONABLE_SURFACE = {"CD33", "CD123", "CD34", "CD117", "CD135", "CD47", "CD38",
                      "CD7", "CD56", "CD25", "CD64"}
SURFACE_KB = {"CD33", "CD123"}     # surface markers with a curated therapy row (flow-pending)


def _num(v):
    try:
        f = float(v)
        return f if f == f else None      # drop NaN
    except (TypeError, ValueError):
        return None


def _truthy(v):
    return str(v).strip().lower() in ("true", "1", "1.0", "yes", "t")


def _clean_marker(c):
    return str(c).split("_")[0].split(".")[0]


def _looks_unknown(name):
    n = str(name).lower()
    return n.startswith("unknown") or n.startswith("unnamed") or n.startswith("nan")


# ---------------------------------------------------------------- LSC witness
def _patient_lsc(ctx, sample_key):
    base = {"witness": "LSC", "kind": "classifier_call"}
    lsc = ctx.tables.get("lsc_calls")
    if lsc is None or sample_key not in lsc.index:
        return {**base, "status": "skipped", "reason": "no LSC call for this sample"}
    row = lsc.loc[sample_key]
    pred = str(row.get("PredictedClass", "")).strip() or None
    maxp = _num(row.get("MaxProb"))
    low, few = _truthy(row.get("LowConfidence")), _truthy(row.get("FewCells"))
    total = _num(row.get("TotalCells"))
    probs = {c: round(_num(row.get(c)) or 0.0, 3) for c in sorted(lsc.columns) if c.startswith("Prob_")}
    confident = bool(pred) and not low and not few
    desc = [f"LSC architecture call: {pred} (MaxProb {round(maxp,3) if maxp is not None else 'NA'}, "
            f"low_confidence={low}, few_cells={few}); upstream RF honest balanced-acc ~0.59 — "
            f"corroborating stemness/risk context, not ground truth."]
    return {**base, "status": "ok", "predicted_class": pred,
            "max_prob": round(maxp, 3) if maxp is not None else None, "class_probs": probs,
            "low_confidence": low, "few_cells": few, "total_cells": total,
            "confident_call": confident,
            "validation_claims": (["lsc"] if confident else []),
            "descriptive_context": desc}


def assess_patient_lsc(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0, "summary": "no LSC call", "caveats": "skipped"}
    ceil = 0.0 if (ev.get("low_confidence") or ev.get("few_cells")) else 0.5
    return _assess_descriptive(client, ev,
        f"You are the LSC (leukemic stem cell architecture) witness. The upstream RF (honest "
        f"balanced accuracy ~0.59, 3-class p-/m-/p+m-LSC) called class {ev.get('predicted_class')} "
        f"at MaxProb {ev.get('max_prob')} (low_confidence={ev.get('low_confidence')}, "
        f"few_cells={ev.get('few_cells')}, probs={ev.get('class_probs')}). This is an RNA-derived "
        f"stemness/risk PHENOTYPE, NOT a genetic subtype.", ceil)


# ---------------------------------------------------------------- surfaceome / ADT witness
def _patient_surfaceome(ctx, sample_key):
    base = {"witness": "surfaceome/ADT", "kind": "imputed_surface"}
    M = dataio.cohort_modality_matrix(ctx, "ADT")
    if M is None or M.shape[0] < 3 or sample_key not in M.index:
        return {**base, "status": "skipped", "reason": "no ADT data / too few samples for a baseline"}
    sd = M.std(ddof=0).replace(0, np.nan)
    z = ((M.loc[sample_key] - M.mean()) / sd).dropna()
    by_marker = {}
    for col, val in z.items():
        by_marker.setdefault(_clean_marker(col), []).append(float(val))
    zmax = {k: round(max(v), 3) for k, v in by_marker.items()}
    high = sorted([(k, v) for k, v in zmax.items() if v >= 1.0], key=lambda kv: (-kv[1], kv[0]))[:10]
    high_markers = [{"marker": k, "z": v} for k, v in high]
    actionable = [m for m in high_markers if m["marker"] in ACTIONABLE_SURFACE]
    therapy_biomarkers = sorted({m["marker"] for m in high_markers if m["marker"] in SURFACE_KB})
    flow_panel = sorted({m["marker"] for m in actionable})
    elevated = ", ".join(f"{m['marker']}(z{m['z']})" for m in high_markers[:6])
    desc = ([f"Elevated (imputed) surface markers vs cohort (z>=1): {elevated}. "
             "Imputed-from-RNA hypotheses; confirm by flow."]
            if high_markers else ["No surface marker is notably elevated (imputed ADT) vs cohort."])
    return {**base, "status": "ok", "high_markers": high_markers, "actionable_markers": actionable,
            "therapy_biomarkers": therapy_biomarkers, "flow_panel": flow_panel,
            "validation_claims": (["surface_marker"] if high_markers else []),
            "descriptive_context": desc}


def assess_patient_surfaceome(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0, "summary": "no ADT data", "caveats": "skipped"}
    return _assess_descriptive(client, ev,
        "You are the surfaceome/ADT (RNA-IMPUTED surface protein) witness. For THIS patient the "
        f"markers elevated vs cohort (z>=1) are {ev.get('high_markers')}; actionable subset "
        f"{ev.get('actionable_markers')}. These are HYPOTHESES requiring flow-cytometry confirmation "
        "before any antigen-directed therapy; ADT is imputed-from-RNA with no held-out fidelity.", 0.4)


# ---------------------------------------------------------------- metabolic / lipid witnesses
def _patient_omics(ctx, sample_key, modality, witness, claim, min_spearman=0.3):
    base = {"witness": witness, "kind": "imputed_omics"}
    M = dataio.cohort_modality_matrix(ctx, modality, min_spearman=min_spearman)
    # need >=3 samples for a robust-z baseline (mirror the surfaceome guard); a 1-2 sample
    # cohort makes MAD all-NaN -> an empty but "ok" read, so skip honestly instead.
    if M is None or M.shape[1] == 0 or M.shape[0] < 3 or sample_key not in M.index:
        return {**base, "status": "skipped", "reason": f"no usable {modality} baseline after fidelity filter"}
    fid = dataio.feature_fidelity(ctx, modality)
    med = M.median()
    mad = (M - med).abs().median().replace(0, np.nan)
    rz = ((M.loc[sample_key] - med) / (1.4826 * mad)).dropna()
    rows, n_named = [], 0
    for feat, zv in rz.items():
        if _looks_unknown(feat):
            continue
        n_named += 1
        f = (_num(fid.get(feat)) if fid is not None else None)
        score = float(zv) * (f if f is not None else 0.0)
        rows.append({"feature": str(feat), "z": round(float(zv), 2),
                     "heldout_spearman": (round(f, 3) if f is not None else None),
                     "weighted_score": round(score, 3)})
    rows.sort(key=lambda r: (-abs(r["weighted_score"]), r["feature"]))
    top = rows[:8]
    fids = [r["heldout_spearman"] for r in top if r["heldout_spearman"] is not None]
    med_fid = round(float(np.median(fids)), 3) if fids else None
    if top:
        desc = [f"{witness}: most-distinctive features (imputed, min held-out Spearman {min_spearman}; "
                f"median fidelity of reported {med_fid}): "
                f"{', '.join(r['feature'] + ('+' if r['z'] > 0 else '-') for r in top[:6])}."]
    else:
        desc = [f"{witness}: no named fidelity-passing feature deviates from cohort for this patient."]
    return {**base, "status": "ok", "min_spearman": min_spearman,
            "n_features_passing_fidelity": int(M.shape[1]),     # cols passing the fidelity gate
            "n_named_features": int(n_named),                   # of those, the named/analyzable ones
            "top_features": top, "median_reported_fidelity": med_fid,
            "validation_claims": ([claim] if top else []),      # don't recommend an assay with no finding
            "descriptive_context": desc}


def _patient_metabolic(ctx, sample_key):
    return _patient_omics(ctx, sample_key, "Metabolite", "metabolic", "metabolite")


def _patient_lipid(ctx, sample_key):
    ev = _patient_omics(ctx, sample_key, "Lipid", "lipid", "lipid_profile")
    if ev.get("status") == "ok":
        cls = {}
        for r in ev["top_features"]:
            c = str(r["feature"]).split("|")[0].split(" ")[0]
            cls[c] = cls.get(c, 0.0) + abs(r["weighted_score"])
        top_cls = sorted(cls.items(), key=lambda kv: (-kv[1], kv[0]))[:4]
        ev["top_lipid_classes"] = [{"lipid_class": c, "score": round(v, 3)} for c, v in top_cls]
    return ev


def assess_patient_metabolic(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0, "summary": "no metabolite data", "caveats": "skipped"}
    return _assess_descriptive(client, ev,
        "You are the metabolic (RNA-IMPUTED metabolomics) witness. For THIS patient the most-distinctive "
        f"fidelity-passing metabolites are {ev.get('top_features')}. Imputed-from-RNA, modest held-out "
        "fidelity (median Spearman ~0.27): descriptive context only, confirm by targeted LC-MS.", 0.35)


def assess_patient_lipid(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0, "summary": "no lipid data", "caveats": "skipped"}
    return _assess_descriptive(client, ev,
        "You are the lipid (RNA-IMPUTED lipidomics) witness. For THIS patient the most-distinctive "
        f"fidelity-passing lipids are {ev.get('top_features')}; enriched lipid classes "
        f"{ev.get('top_lipid_classes')}. Imputed-from-RNA, modest fidelity (median Spearman ~0.36): "
        "descriptive only, confirm by untargeted lipidomics.", 0.35)


# ---------------------------------------------------------------- GRN witness
def _patient_grn(ctx, sample_key):
    base = {"witness": "GRN-regulon", "kind": "imputed_regulon"}
    row = dataio.sample_modality_matrix(ctx, "GRN", [sample_key])
    if row is None or row.empty:
        return {**base, "status": "skipped", "reason": "no GRN data for this sample"}
    r = row.iloc[0]
    a = ctx.open_modality("GRN")
    tfmap = (dict(zip(a.var_names.astype(str), a.var["TF"].astype(str)))
             if "TF" in a.var.columns else {})
    agg = {}
    for col, val in r.items():
        agg[tfmap.get(str(col), str(col).split("|")[0])] = \
            agg.get(tfmap.get(str(col), str(col).split("|")[0]), 0.0) + float(val)
    top_tfs = [{"tf": t, "activity": round(v, 3)}
               for t, v in sorted(agg.items(), key=lambda kv: (-kv[1], kv[0]))[:8]]
    top_edges = [{"edge": str(c), "activity": round(float(v), 3)}
                 for c, v in sorted(r.items(), key=lambda kv: (-float(kv[1]), str(kv[0])))[:8]]
    desc = [f"Top active regulators (imputed GRN — fidelity UNKNOWN, no held-out metric): "
            f"{', '.join(t['tf'] for t in top_tfs[:5])}. Discovery context only."]
    return {**base, "status": "ok", "fidelity_status": "unknown",
            "top_active_tfs": top_tfs, "top_active_edges": top_edges,
            "descriptive_context": desc}


def assess_patient_grn(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0, "summary": "no GRN data", "caveats": "skipped"}
    return _assess_descriptive(client, ev,
        "You are the GRN (gene regulatory network, RNA-IMPUTED TF->target activity) witness. For THIS "
        f"patient the top active regulators are {ev.get('top_active_tfs')}. CRITICAL: GRN has NO held-out "
        "fidelity metric (unknown reliability) — discovery context only, never a subtype or therapy call.", 0.3)


# ---------------------------------------------------------------- cell-communication witness
def _patient_signaling(ctx, sample_key):
    base = {"witness": "cell-communication", "kind": "signaling"}
    row = dataio.cellcomm_matrix(ctx, [sample_key])
    if row is None or row.empty:
        return {**base, "status": "skipped", "reason": "no cell-communication data for this sample"}
    r = row.iloc[0]
    nz = [(str(c), float(v)) for c, v in r.items() if float(v) > 0]
    top_int = sorted(nz, key=lambda kv: (-kv[1], kv[0]))[:12]
    axes, senders, receivers = {}, {}, {}
    for c, v in nz:
        parts = c.split("|")
        if len(parts) == 4:
            s_, lig, rec, rcv = parts
            ax = f"{lig}->{rec}"
            axes[ax] = max(axes.get(ax, 0.0), v)
            senders[s_] = senders.get(s_, 0.0) + v
            receivers[rcv] = receivers.get(rcv, 0.0) + v
    top_axes = [{"axis": a, "score": round(v, 3)}
                for a, v in sorted(axes.items(), key=lambda kv: (-kv[1], kv[0]))[:10]]
    if top_axes:
        desc = [f"Dominant signaling axes (independent fastComm L-R inference): "
                f"{', '.join(a['axis'] for a in top_axes[:5])}; {len(nz)} called interactions."]
    elif nz:   # interaction names not in the expected sender|ligand|receptor|receiver form
        desc = [f"{len(nz)} called interactions, but var names are not in the expected "
                f"4-part sender|ligand|receptor|receiver form; top: "
                f"{', '.join(c for c, _ in top_int[:3])}."]
    else:
        desc = ["no called ligand-receptor interactions for this patient."]
    return {**base, "status": "ok", "n_edges": len(nz),
            "top_interactions": [{"interaction": c, "score": round(v, 3)} for c, v in top_int],
            "top_axes": top_axes,
            "top_sender_states": [s for s, _ in sorted(senders.items(), key=lambda kv: (-kv[1], kv[0]))[:6]],
            "top_receiver_states": [s for s, _ in sorted(receivers.items(), key=lambda kv: (-kv[1], kv[0]))[:6]],
            "descriptive_context": desc}


def assess_patient_signaling(client, ev):
    if ev.get("status") != "ok":
        return {"confidence": 0.0, "reliability_weight": 0.0, "summary": "no cell-comm data", "caveats": "skipped"}
    return _assess_descriptive(client, ev,
        "You are the cell-communication witness (fastComm ligand-receptor inference — a relatively "
        f"INDEPENDENT axis, derived from measured cell-frequency + expression). For THIS patient the "
        f"dominant signaling axes are {ev.get('top_axes')} over {ev.get('n_edges')} called interactions; "
        f"hub senders {ev.get('top_sender_states')}, receivers {ev.get('top_receiver_states')}. "
        "Descriptive signaling context, NOT a subtype call.", 0.5)


def _assess_descriptive(client, ev, role, weight_ceiling):
    """Shared LLM assessment for the descriptive witnesses: narrate + self-cap the weight.
    The weight is reporting-only (these domains never vote), but kept honest."""
    ctxstr = "; ".join(ev.get("descriptive_context", []))[:700]
    if client is not None:
        try:
            out = client.chat_json(
                f"{role}\nSummarize concisely what this means FOR THIS PATIENT and how much to trust it. "
                f"This is CORROBORATING/DESCRIPTIVE evidence (reliability_weight must be <= {weight_ceiling}); "
                "it must not be treated as a subtype call or as independent confirmation of an RNA-based read. "
                'Return ONLY JSON {"confidence":0..1,"reliability_weight":0..1,"summary":str,"caveats":str}.',
                required=OPINION_KEYS)
            out["reliability_weight"] = round(min(float(out.get("reliability_weight") or 0.0),
                                                  weight_ceiling), 2)
            return out
        except (LLMError, Exception):
            pass
    return {"confidence": round(min(0.5, weight_ceiling + 0.1), 2),
            "reliability_weight": round(weight_ceiling, 2),
            "summary": (ctxstr[:200] or role[:120]),
            "caveats": "LLM unavailable; deterministic fallback (descriptive, capped weight)"}


def reconcile_patient(client, sample_key, items):
    brief = [{"witness": it["evidence"].get("witness"), "kind": it["evidence"].get("kind"),
              "confidence": it["opinion"].get("confidence"),
              "reliability_weight": it["opinion"].get("reliability_weight"),
              "summary": it["opinion"].get("summary")} for it in items]
    prompt = ("You are the chair producing a PER-PATIENT AML triage report. Reconcile these witness reads "
              f"for one patient into a decision:\n{json.dumps(brief, default=str)[:3800]}\n"
              "Weight INDEPENDENT axes (the genetic witness, measured composition) and held-out / above-chance "
              "predictive evidence; treat imputed-from-RNA as corroborating, NOT independent. Give the "
              "best-supported subtype call (with confidence), the targetable therapies worth considering FOR "
              "THIS PATIENT (with rationale), concrete recommended validations (flow markers, targeted "
              "sequencing, ex-vivo drug test), any conflicts, and a short rationale. Be candid about uncertainty. "
              'Return ONLY JSON {"subtype_call":str,"targetable_therapies":str,"recommended_validations":str,'
              '"overall_confidence":"low"|"medium"|"high","conflicts":str,"rationale":str}.')
    try:
        return client.chat_json(prompt, required=PATIENT_CONSENSUS_KEYS)
    except (LLMError, Exception) as e:
        return {"subtype_call": "(LLM unavailable)", "targetable_therapies": "",
                "recommended_validations": "", "overall_confidence": "low", "conflicts": "",
                "rationale": str(e)[:200]}


def run_patient_panel(ctx, sample_key, blocks=None, with_genetic=True, with_udon=True,
                      strategy="donor_kfold", permutations=40, client=None,
                      max_rounds=0, mode="continuous"):
    """Drive the agent roster → shared ledger → arbiter v2 (deterministic anchor) for one patient.
    max_rounds=0 -> single pass (default). max_rounds>=1 -> Phase C guarded feedback loop: witnesses
    deterministically defer toward the observed genetic anchor over bounded rounds (the anchor still
    has final say; round-0 baseline + groupthink drift recorded). mode = continuous | conflict_triggered."""
    from . import ledger as _led, arbiter
    from .agents import default_roster
    from .agent import run_agent
    client = client or LLMClient()
    scope = {"mode": "patient", "sample_key": sample_key,
             "strategy": strategy, "permutations": permutations}
    led = _led.Ledger(sample_key, "subtype", provenance=ctx.getart("provenance"))
    ctx.ledger = led
    for spec in default_roster(blocks, with_genetic, with_udon):
        res = run_agent(spec, ctx, scope, client)
        led.append(res, round=0)
    led.persist(ctx)
    consensus = arbiter.reconcile_patient(client, ctx, led)
    led.set_arbiter(consensus, round=0)
    deliberation = None
    if max_rounds and int(max_rounds) >= 1:
        from . import feedback
        deliberation = feedback.deliberate(ctx, led, client, max_rounds=int(max_rounds), mode=mode)
        consensus = deliberation["consensus"]
    else:
        led.finalize("single_pass")
    led.persist(ctx)
    s = ctx.tables["samples"]
    srow = s.loc[sample_key] if sample_key in s.index else None
    report = {"mode": "patient_panel", "sample_key": sample_key,
              "annotation": (str(srow["annotation"]) if srow is not None else None),
              "dataset": (str(srow["dataset"]) if srow is not None else None),
              "provenance": ctx.getart("provenance"),
              "panel": [{"witness": e["witness"], "grounding": e["grounding"],
                         "independence": e["independence"], "evidence": e["evidence"],
                         "opinion": e["opinion"]} for e in led.current_entries()],
              "consensus": consensus, "deliberation": deliberation}
    ctx.save_json(report, "patient_report.json")
    _write_patient_md(ctx, report)
    return report


def _write_patient_md(ctx, report):
    c = report["consensus"]
    L = [f"# Patient triage — {report['sample_key']}",
         f"(atlas annotation: {report.get('annotation')} · dataset: {report.get('dataset')} · "
         f"KB {c.get('knowledge_version')})", "", "## Witness reads"]
    for it in report["panel"]:
        o = it["opinion"]
        L.append(f"- **{it.get('witness')}** ({it.get('independence')} / {it.get('grounding')}, "
                 f"conf {o.get('confidence')}, weight {o.get('reliability_weight')}): {o.get('summary')}")
    L += ["", f"## Decision ({c.get('overall_confidence')} confidence)",
          f"- **Subtype call:** {c.get('subtype_call')}  (concordance {c.get('concordance')})",
          f"- **Per-witness consistency:** {c.get('per_witness_consistency')}",
          f"- **Conflicts:** {c.get('conflicts') or 'none'}", "",
          "### Ranked therapy hypotheses (knowledge-grounded)"]
    for t in (c.get("ranked_therapy_hypotheses") or []):
        L.append(f"- {t.get('biomarker')} → {t.get('drug')}  [{t.get('evidence_level')}; {t.get('source')}]")
    if not c.get("ranked_therapy_hypotheses"):
        L.append("- (none — no observed targetable driver)")
    surf = c.get("surface_therapy_hypotheses") or []
    if surf:
        L += ["", "### Surface-marker therapy hypotheses (imputed ADT — FLOW-PENDING, not orders)"]
        for t in surf:
            L.append(f"- {t.get('biomarker')} → {t.get('drug')}  "
                     f"[{t.get('evidence_level')}; {t.get('source')}; requires flow confirmation]")
    L += ["", "### Recommended validations"]
    for v in (c.get("recommended_validations") or []):
        L.append(f"- {v.get('claim')}: {v.get('validation')}")
    desc = c.get("descriptive_findings") or []
    if desc:
        L += ["", "### Descriptive / discovery context (corroborating; non-voting)"]
        for d in desc:
            L.append(f"- _{d.get('witness')}_ ({d.get('domain')}): {d.get('note')}")
    dl = report.get("deliberation")
    if dl:
        d = dl.get("drift", {})
        L += ["", f"### Deliberation (Phase C, {dl.get('mode')})",
              f"- **Deliberation rounds:** {dl.get('deliberation_rounds')} ({dl.get('stop_reason')})",
              f"- **Baseline → final leading:** {d.get('baseline_leading')} → {d.get('final_leading')} "
              f"(changed: {d.get('leading_changed')}; final genetically confirmed: "
              f"{d.get('final_genetically_confirmed')})",
              f"- **Concordance baseline → final:** {d.get('baseline_concordance')} → "
              f"{d.get('final_concordance')}",
              f"- **Groupthink warning:** {d.get('groupthink_warning')}"]
    L += ["", "### Rationale", c.get("rationale", "")]
    fp = ctx.path("PATIENT.md")
    with open(fp, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return fp


def _write_md(ctx, report):
    L = [f"# MOSAIC-AML panel — {report['target']}", "",
         "| witness | kind | bal.acc | perm p | confidence | weight |",
         "|---|---|---|---|---|---|"]
    for it in report["panel"]:
        e, o = it["evidence"], it["opinion"]
        L.append(f"| {e.get('witness')} | {e.get('kind')} | {e.get('balanced_accuracy','—')} | "
                 f"{e.get('permutation_pvalue','—')} | {o.get('confidence')} | {o.get('reliability_weight')} |")
    c = report["consensus"]
    L += ["", f"**Consensus ({c.get('overall_confidence')} confidence):** {c.get('consensus')}", "",
          f"**Targetable summary:** {c.get('targetable_summary')}", "",
          f"**Recommended validations:** {c.get('recommended_validations')}", "",
          f"**Per-witness consistency:** {c.get('per_witness_consistency')}", "",
          f"**Conflicts:** {c.get('conflicts')}", "",
          f"_Caveat:_ {c.get('caveat')}"]
    fp = ctx.path("PANEL.md")
    with open(fp, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return fp
