"""Therapy hypotheses + confirmatory tests — closes the loop from a MOSAIC-AML call to an action.

RESEARCH DECISION SUPPORT ONLY. Everything here is a hypothesis generated from drivers that were
either PREDICTED from expression (never sequenced — must be confirmed before they mean anything) or
SUPPLIED as known calls. This is not medical advice and is not a prescribing tool.

Two panels are produced from a report's mutation calls:
  * treatments — per driver, the targeted agents that driver implicates, with an evidence tier and a
    clear predicted-vs-confirmed status. Keyed at the VARIANT level where the variant changes the drug
    (e.g. FLT3-ITD vs FLT3-TKD: type-II inhibitors like quizartinib do not cover TKD; gilteritinib does).
  * tests — what to order to confirm/refute the calls, prioritised. Deliberately also recommends testing
    HIGH-VALUE actionable drivers the model is UNRELIABLE about (abstain / low AUROC): the point of the
    loop is to sequence exactly what the model can't settle.

The loop closes by feeding confirmed results back in via `--mutations`, which engages the deterministic
genetic anchor on the next run.
"""
from __future__ import annotations
import json
import os

_D = lambda name, cls, tier, note="": {"drug": name, "class": cls, "tier": tier, "note": note}

# ---- variant-level rules (the variant changes the answer) -----------------------------------------
BY_CATEGORY = {
    "FLT3_ITD": {
        "drugs": [_D("gilteritinib", "FLT3 inhibitor (type I)", "approved · R/R", "covers ITD and TKD"),
                  _D("midostaurin", "FLT3 inhibitor (type I)", "approved · newly dx + chemo"),
                  _D("quizartinib", "FLT3 inhibitor (type II)", "approved · newly dx, ITD", "ITD only — no TKD activity")],
        "test": "FLT3-ITD PCR with allelic ratio",
        "why": "ITD allelic ratio guides inhibitor choice and transplant decisions",
    },
    "FLT3_TKD_D835/I836": {
        "drugs": [_D("gilteritinib", "FLT3 inhibitor (type I)", "approved · R/R", "type I — active against TKD"),
                  _D("midostaurin", "FLT3 inhibitor (type I)", "approved · newly dx + chemo")],
        "test": "FLT3 TKD (D835/I836) sequencing",
        "why": "type-II inhibitors (quizartinib) are NOT active against TKD — the variant changes the drug",
    },
    "FLT3_other_TKD_or_JM": {
        "drugs": [_D("gilteritinib", "FLT3 inhibitor (type I)", "approved · R/R")],
        "test": "FLT3 full-gene NGS", "why": "non-canonical FLT3 activation; confirm the exact residue",
    },
    "IDH1_R132": {
        "drugs": [_D("ivosidenib", "IDH1 inhibitor", "approved"),
                  _D("olutasidenib", "IDH1 inhibitor", "approved")],
        "test": "IDH1 R132 hotspot NGS", "why": "IDH1 inhibitors require a confirmed R132 variant",
    },
    "IDH2_R140": {
        "drugs": [_D("enasidenib", "IDH2 inhibitor", "approved")],
        "test": "IDH2 R140/R172 hotspot NGS", "why": "IDH2 inhibitor requires a confirmed variant",
    },
    "IDH2_R172": {
        "drugs": [_D("enasidenib", "IDH2 inhibitor", "approved")],
        "test": "IDH2 R140/R172 hotspot NGS", "why": "IDH2 inhibitor requires a confirmed variant",
    },
    "NPM1_exon12_frameshift": {
        "drugs": [_D("revumenib / ziftomenib", "menin inhibitor", "approved (R/R) · NPM1-mutant")],
        "test": "NPM1 exon-12 PCR", "why": "also the best MRD marker in NPM1-mutant AML",
    },
    "KMT2A_fusion": {
        "drugs": [_D("revumenib", "menin inhibitor", "approved · R/R KMT2A-rearranged")],
        "test": "karyotype + KMT2A FISH / fusion RNA panel", "why": "identifies the partner and confirms rearrangement",
    },
    "KMT2A_PTD": {
        "drugs": [_D("menin inhibitor", "menin inhibitor", "investigational · KMT2A-PTD")],
        "test": "KMT2A-PTD assay", "why": "PTD is missed by standard SNV panels",
    },
    "KIT_D816": {
        "drugs": [_D("avapritinib", "KIT inhibitor", "approved (SM) · CBF-AML investigational"),
                  _D("dasatinib", "multi-kinase / KIT", "investigational · CBF-AML")],
        "test": "KIT D816 NGS (CBF-AML context)", "why": "KIT mutation worsens CBF-AML prognosis",
    },
    "JAK2_V617F": {
        "drugs": [_D("ruxolitinib", "JAK1/2 inhibitor", "approved (MPN) · secondary-AML context")],
        "test": "JAK2 V617F PCR", "why": "suggests an MPN-derived / secondary AML",
    },
    "TP53_hotspot_DBD": {
        "drugs": [_D("HMA + venetoclax", "hypomethylating + BCL2", "standard of care · poor response"),
                  _D("clinical trial", "—", "preferred", "no approved TP53-directed agent")],
        "test": "TP53 NGS + karyotype/FISH (allelic state)", "why": "adverse risk; allelic state (mono vs bi) drives prognosis",
    },
    "TP53_LOF/splice/frameshift": {
        "drugs": [_D("HMA + venetoclax", "hypomethylating + BCL2", "standard of care · poor response"),
                  _D("clinical trial", "—", "preferred", "no approved TP53-directed agent")],
        "test": "TP53 NGS + karyotype/FISH (allelic state)", "why": "adverse risk; commonly with complex karyotype",
    },
}

# splicing factors — the drug-repositioning lead (SpliceScout convergence)
for _c in ("SF3B1_K700", "SF3B1_K666", "U2AF1_S34", "U2AF1_Q157/R156", "SRSF2", "ZRSR2"):
    BY_CATEGORY[_c] = {
        "drugs": [_D("splicing modulator (e.g. H3B-8800 class)", "spliceosome modulator", "investigational",
                     "drug-repositioning lead — spliceosome-mutant cells are selectively sensitive")],
        "test": "splicing-factor NGS (SF3B1 / SRSF2 / U2AF1 / ZRSR2)",
        "why": "spliceosome-mutant AML is the repositioning target — pair with a splicing readout",
    }

# ---- gene-level fallback (variant doesn't change the answer) --------------------------------------
_RAS = {
    "drugs": [_D("trametinib", "MEK inhibitor", "investigational · RAS-pathway"),
              _D("SHP2 inhibitor", "SHP2", "investigational")],
    "test": "RAS-pathway NGS (NRAS/KRAS/PTPN11/NF1/CBL)",
    "why": "RAS-pathway activation is linked to venetoclax resistance and relapse",
}
BY_GENE = {
    "NRAS": _RAS, "KRAS": _RAS, "PTPN11": _RAS, "NF1": _RAS, "CBL": _RAS,
    "DNMT3A": {"drugs": [_D("hypomethylating agent", "HMA", "context / investigational")],
               "test": "DNMT3A NGS (R882 vs non-R882)",
               "why": "R882 is the functionally distinct hotspot — worth resolving"},
    "ASXL1": {"drugs": [_D("standard of care", "—", "prognostic")], "test": "myeloid NGS panel",
              "why": "adverse-risk marker (ELN); no targeted agent"},
    "RUNX1": {"drugs": [_D("standard of care", "—", "prognostic")], "test": "myeloid NGS panel",
              "why": "adverse-risk marker (ELN); no targeted agent"},
    "TET2": {"drugs": [_D("hypomethylating agent", "HMA", "context / investigational")],
             "test": "myeloid NGS panel", "why": "clonal-haematopoiesis driver; context for HMA"},
    "WT1": {"drugs": [_D("standard of care", "—", "prognostic")], "test": "myeloid NGS panel + WT1 MRD",
            "why": "WT1 transcript can serve as an MRD marker"},
}

# always-order baseline workup
BASELINE_TESTS = [
    {"test": "Karyotype + FISH", "why": "cytogenetic risk (inv16, t(8;21), -5/-7, complex) — bulk RNA does not call these",
     "priority": "baseline"},
    {"test": "Targeted myeloid NGS panel", "why": "the ground truth this panel's predictions must be confirmed against",
     "priority": "baseline"},
    {"test": "Flow cytometry immunophenotype", "why": "blast %, lineage, and an MRD baseline", "priority": "baseline"},
]

ABSTAIN_AUC = 0.65        # below this the model is unreliable for that category


def _gene_of(cat):
    return str(cat).split("_")[0].split("-")[0].upper()


def _rule_for(cat):
    r = BY_CATEGORY.get(cat)
    if r:
        return r, "variant"
    r = BY_GENE.get(_gene_of(cat))
    return (r, "gene") if r else (None, None)


def _fallback_explanation(t, co_drivers):
    """Deterministic explanation — always present, and what you get with --no-llm."""
    drugs = ", ".join(d["drug"] for d in t.get("drugs", []))
    conf = t.get("confidence")
    auc = t.get("reliability_auc")
    if t.get("predicted"):
        basis = ("This driver was PREDICTED from expression (model score %s, reliability AUROC %s) and has "
                 "NOT been sequenced — %s is implicated only if the call is confirmed."
                 % ("%.2f" % conf if isinstance(conf, (int, float)) else "n/a",
                    "%.2f" % auc if isinstance(auc, (int, float)) else "n/a", drugs))
    else:
        basis = "This driver was sequenced and confirmed, which is what implicates %s." % drugs
    why = t.get("rationale") or ""
    co = ""
    if co_drivers:
        co = " Co-occurring calls on this sample: %s." % ", ".join(co_drivers[:4])
    return (basis + (" " + why if why else "") + co).strip()


def _agent_json(client, prompt, max_tokens):
    """Ask for JSON WITHOUT response_format, then parse leniently.

    llm.py's chat_json() forces response_format={"type":"json_object"}; vllm/gpt-oss stalls on that
    guided decoding and the request dies at the 180s timeout. Plain completions come back in ~2s and do
    emit valid JSON when asked for it — so ask plainly and tolerate fences / trailing prose.
    """
    txt = client.chat(prompt, max_tokens=max_tokens)
    s = (txt or "").strip()
    if "```" in s:                                    # strip a ```json fence if present
        parts = s.split("```")
        for p in parts:
            p = p.strip()
            if p.lower().startswith("json"):
                p = p[4:].strip()
            if p.startswith("{"):
                s = p
                break
    i, j = s.find("{"), s.rfind("}")                  # tolerate prose around the object
    if i >= 0 and j > i:
        s = s[i:j + 1]
    return json.loads(s)


def explain_treatments(client, treatments, context=None):
    """Agent narration for WHY each implicated therapy follows from the driver.

    Constrained exactly like the arbiter's narrator: the deterministic rules ALREADY chose the drugs.
    The agent may only EXPLAIN them — it must not add, drop or substitute an agent, must not invent
    evidence, and must not treat a PREDICTED driver as confirmed. Every treatment gets a deterministic
    explanation first, so a missing/failed LLM degrades to a real answer rather than nothing.
    """
    if not treatments:
        return treatments
    co = [t["driver"] for t in treatments]
    for t in treatments:                                   # deterministic baseline for every row
        t["explanation"] = _fallback_explanation(t, [c for c in co if c != t["driver"]])
        t["explained_by"] = "rules"
    if client is None:
        return treatments

    brief = [{"driver": t["driver"], "status": t["status"], "model_score": t.get("confidence"),
              "reliability_auroc": t.get("reliability_auc"),
              "drugs": [d["drug"] for d in t.get("drugs", [])],
              "drug_tiers": [d.get("tier") for d in t.get("drugs", [])],
              "rule_rationale": t.get("rationale")} for t in treatments]
    allowed = {t["driver"]: {d["drug"] for d in t.get("drugs", [])} for t in treatments}
    try:
        out = _agent_json(client,
            "You are the AML tumor-board pharmacologist. A deterministic rule engine has ALREADY decided "
            "which agents each driver implicates; your ONLY job is to explain WHY, for a clinician. "
            "Panel (with this patient's evidence): %s. %s"
            "For each driver write 2-3 sentences: what the driver does biologically, why THOSE agents "
            "follow from it, and how much weight to put on it given the evidence shown. "
            "HARD RULES: (1) never name an agent that is not in that driver's drug list; (2) never add, "
            "drop or substitute a therapy; (3) if status says PREDICTED, state plainly that it is "
            "inferred from expression and must be confirmed by sequencing before it guides treatment; "
            "(4) invent no evidence, numbers or trial claims beyond what is given. "
            'Return ONLY JSON {"explanations": {"<driver>": "<text>"}}.'
            % (json.dumps(brief, default=str)[:3000],
               ("Sample context: %s. " % json.dumps(context, default=str)[:400]) if context else ""),
            # bounded budget: llm.py defaults to 500k, which vllm will not finish
            int(os.environ.get("AMLMM_THERAPY_MAX_TOKENS", "3000")))
        got = out.get("explanations") or {}
        for t in treatments:
            txt = got.get(t["driver"])
            if not isinstance(txt, str) or not txt.strip():
                continue
            # honesty guard: reject narration that names a drug outside this driver's decided list
            low = txt.lower()
            other = set()
            for d2, ds in allowed.items():
                if d2 == t["driver"]:
                    continue
                other |= ds
            leaked = [x for x in other - allowed[t["driver"]] if x.lower() in low and len(x) > 4]
            if leaked:
                continue                                   # keep the deterministic text instead
            if t.get("predicted") and "confirm" not in low and "sequenc" not in low:
                txt = txt.strip() + " (Predicted from expression — confirm by sequencing before acting.)"
            t["explanation"] = txt.strip()[:900]
            t["explained_by"] = "agent"
    except Exception:
        pass                                               # deterministic explanations already in place
    return treatments


def build_panels(predictions, supplied=None, client=None, context=None):
    """predictions: the report's mutation_predictions list (each: mutation, call, probability,
    confidence, heldout_auc/cv_auroc). supplied: user-provided (sequenced) drivers, gene-level.
    -> {"treatments": [...], "tests": [...], "note": ...}"""
    preds = list(predictions or [])
    supplied = [str(s).upper() for s in (supplied or [])]
    sup_set = set(supplied)

    def _exact(cat):                          # the exact variant category was supplied
        return str(cat).upper() in sup_set

    def _gene_sup(cat):                       # the GENE was supplied, but not which variant class
        return _gene_of(cat) in sup_set

    def _confirmed(cat):                      # this driver's gene was actually sequenced
        return _exact(cat) or _gene_sup(cat)

    treatments, tests, seen_t, seen_x = [], [], set(), set()

    # 1) drivers the model CALLS PRESENT (confident), or whose exact variant was SUPPLIED -> therapy
    #    NB a GENE-level supply (e.g. "TP53") must NOT conjure a row for a variant class the model
    #    never called - sequencing the gene tells you it's mutated, not that every class of it is.
    for p in preds:
        cat = p.get("mutation")
        conf_ok = str(p.get("confidence", "")) == "ok"
        called = p.get("call") == "present"
        if not ((called and conf_ok) or _exact(cat)):
            continue
        rule, level = _rule_for(cat)
        if not rule:
            continue
        key = (cat,)
        if key in seen_t:
            continue
        seen_t.add(key)
        if _exact(cat):
            status = "confirmed (sequenced)"
        elif _gene_sup(cat):
            status = "gene confirmed · variant predicted"
        else:
            status = "PREDICTED — confirm before acting"
        treatments.append({
            "driver": cat, "match": level, "status": status,
            "predicted": not _exact(cat),
            "confidence": p.get("probability"), "reliability_auc": p.get("heldout_auc", p.get("cv_auroc")),
            "drugs": rule["drugs"], "rationale": rule.get("why", ""),
        })
        # the matching confirmatory test, highest priority
        t = rule.get("test")
        if t and t not in seen_x:
            seen_x.add(t)
            tests.append({"test": t, "why": rule.get("why", ""), "driver": cat,
                          "priority": "confirm" if not _confirmed(cat) else "done",
                          "reason": ("confirm the predicted driver before any targeted decision"
                                     if not _confirmed(cat) else "already sequenced — no action")})

    # 2) actionable drivers the model is UNRELIABLE about -> sequence them anyway (this is the point)
    for p in preds:
        cat = p.get("mutation")
        if _confirmed(cat):
            continue
        auc = p.get("heldout_auc", p.get("cv_auroc"))
        weak = (auc is None) or (isinstance(auc, (int, float)) and auc < ABSTAIN_AUC) \
            or str(p.get("confidence", "")).startswith("abstain")
        if not weak:
            continue
        rule, _ = _rule_for(cat)
        if not rule or cat not in BY_CATEGORY:      # only chase genuinely actionable variants
            continue
        t = rule.get("test")
        if not t or t in seen_x:
            continue
        seen_x.add(t)
        tests.append({"test": t, "why": rule.get("why", ""), "driver": cat, "priority": "resolve",
                      "reason": "actionable driver the model cannot settle (AUROC %s) — sequence it"
                                % ("n/a" if auc is None else round(float(auc), 2))})

    # 3) baseline workup
    for b in BASELINE_TESTS:
        if b["test"] not in seen_x:
            seen_x.add(b["test"])
            tests.append({**b, "driver": None, "reason": "standard workup"})

    order = {"confirm": 0, "resolve": 1, "baseline": 2, "done": 3}
    tests.sort(key=lambda x: order.get(x.get("priority"), 9))
    treatments.sort(key=lambda t: (t["status"].startswith("PREDICTED"), -(t.get("confidence") or 0)))
    explain_treatments(client, treatments, context)        # agent explains; rules already decided

    return {
        "treatments": treatments,
        "tests": tests,
        "note": "RESEARCH DECISION SUPPORT — not medical advice. Drivers marked PREDICTED were inferred "
                "from expression and have NOT been sequenced; confirm them before any targeted decision. "
                "Close the loop: re-run with --mutations <confirmed drivers> to engage the deterministic "
                "genetic anchor.",
    }
