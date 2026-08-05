"""Phase D — retrospective clinical validation (honest, underpowered-aware).

Question: does the engine's REASONING track independent clinical reality? We take the
engine's DETERMINISTIC cohort-wide outputs — the genetic-anchored driver (the arbiter's
leading hypothesis for anchored patients) and its ELN-expected risk class, plus the upstream
LSC stemness call — and test their ASSOCIATION with the sparse clinical labels actually
present (ELN_risk, clinical_response, overall_survival).

These are RETROSPECTIVE ASSOCIATIONS on tiny, partially-overlapping label sets — NOT trained
predictors, NOT validated models. Every result carries its n and an `underpowered` flag
(n < 20); a null association is a legitimate, reportable finding. Honesty over flattery.

Stats: 2x2 -> Fisher exact (odds ratio + exact p); k-class concordance -> agreement rate +
a seeded label-permutation p (matches the pipeline's permutation-baseline ethos). No asymptotic
chi-square on tiny cells.
"""
from __future__ import annotations
import numpy as np

from . import genetics
from .arbiter import ANCHOR_MAP, ANCHOR_PRIORITY

UNDERPOWERED_N = 20

# Heuristic ELN-2022-flavored expected risk for the SINGLE anchored driver. A simplification:
# real ELN integrates co-mutations (e.g. NPM1+FLT3-ITD is intermediate) — but the anchor's
# priority already routes NPM1+FLT3 to FLT3, so the single-driver proxy is a defensible test.
ELN_EXPECTED = {
    # favorable
    "NPM1": "Favorable", "Inv16": "Favorable", "t(8;21)": "Favorable", "CEBPA": "Favorable",
    # adverse
    "TP53": "Adverse", "Complex": "Adverse", "DEL(7)": "Adverse", "DEL(5)": "Adverse",
    "KMT2Ar": "Adverse", "ASXL1": "Adverse", "RUNX1": "Adverse", "SF3B1": "Adverse",
    "SRSF2": "Adverse", "U2AF1": "Adverse", "APL": "Favorable",
    # intermediate / not risk-defining alone
    "FLT3": "Intermediate", "IDH1": "Intermediate", "IDH2": "Intermediate",
    "DNMT3A": "Intermediate", "TET2": "Intermediate", "NRAS": "Intermediate",
    "KRAS": "Intermediate", "WT1": "Intermediate", "CKIT": "Intermediate", "KIT": "Intermediate",
    "CSF3R": "Intermediate", "Trisomy8": "Intermediate",
}


def _anchored_driver(M, row):
    present = [c.replace("mut_", "").replace("cyto_", "")
               for c in M.filter(regex="^(mut_|cyto_)").columns if row.get(c) == 1.0]
    for d in sorted(present, key=lambda d: ANCHOR_PRIORITY.index(d) if d in ANCHOR_PRIORITY else 99):
        if d in ANCHOR_MAP:
            return ANCHOR_MAP[d]
    return None


def _norm_eln(v):
    if v is None:
        return None
    t = str(v).strip().lower()
    if t in ("", "nan", "unknown", "na"):
        return None
    if "favor" in t:
        return "Favorable"
    if "adverse" in t or "poor" in t or "high" in t:
        return "Adverse"
    if "interm" in t:
        return "Intermediate"
    return None


def _responder(v):
    """Binarize the (messy) clinical_response column. responder=True if it carries a
    '(responder)' tag; False if it reads 'adverse'/'refractory'/'non'; else None (excluded)."""
    if v is None:
        return None
    t = str(v).strip().lower()
    if t in ("", "nan", "unknown", "na"):
        return None
    if "responder" in t and "non" not in t:
        return True
    if "adverse" in t or "refractor" in t or "resistant" in t or t.startswith("non"):
        return False
    return None


def build_table(ctx):
    import pandas as pd
    s = ctx.tables["samples"]
    M = ctx.tables.get("mutations")
    if M is None:
        M = genetics.build_mutation_matrix(ctx)
    anchored = M.apply(lambda r: _anchored_driver(M, r), axis=1)
    lsc = ctx.tables.get("lsc_calls")
    t = pd.DataFrame(index=s.index)
    t["anchored_driver"] = anchored.reindex(s.index)
    t["expected_eln"] = t["anchored_driver"].map(lambda d: ELN_EXPECTED.get(d) if d else None)
    t["eln_actual"] = s["ELN_risk"].map(_norm_eln) if "ELN_risk" in s.columns else None
    t["responder"] = s["clinical_response"].map(_responder) if "clinical_response" in s.columns else None
    t["survival"] = pd.to_numeric(s["overall_survival"], errors="coerce") if "overall_survival" in s.columns else np.nan
    if lsc is not None and "PredictedClass" in lsc.columns:
        low = lsc["LowConfidence"].astype(str).str.lower().isin(["true", "1", "1.0", "yes"]) if "LowConfidence" in lsc.columns else False
        cls = lsc["PredictedClass"].where(~low)
        t["lsc_class"] = cls.reindex(s.index)
    else:
        t["lsc_class"] = None
    return t


def _fisher2x2(bool_a, bool_b):
    """Fisher exact on two boolean Series. Drops any null rows first (defensive); callers
    should pass already-clean booleans on an explicit support."""
    from scipy.stats import fisher_exact
    import pandas as pd
    df = pd.DataFrame({"a": bool_a, "b": bool_b}).dropna()
    if df.empty:
        return {"n": 0, "table": [[0, 0], [0, 0]], "odds_ratio": None, "p": None, "underpowered": True}
    aa, bb = df["a"].astype(bool), df["b"].astype(bool)
    a = int((aa & bb).sum()); b = int((aa & ~bb).sum())
    c = int((~aa & bb).sum()); d = int((~aa & ~bb).sum())
    try:
        orr, p = fisher_exact([[a, b], [c, d]])
    except Exception:
        orr, p = None, None
    n = a + b + c + d
    return {"n": n, "table": [[a, b], [c, d]],
            "odds_ratio": (round(float(orr), 3) if orr is not None and np.isfinite(orr) else None),
            "p": (round(float(p), 5) if p is not None else None), "underpowered": n < UNDERPOWERED_N}


def _assoc(t, col_a, pred_a, col_b, pred_b):
    """2x2 association on the EXPLICIT both-columns-present support. `.notna()` handles None
    AND NaN consistently across pandas versions (a pre-built mask + `x is not None` does NOT
    — it leaked non-anchored rows in as False on the cluster), so n / odds ratios are
    environment-independent and reproducible."""
    sub = t[t[col_a].notna() & t[col_b].notna()]
    if len(sub) == 0:
        return {"n": 0, "table": [[0, 0], [0, 0]], "odds_ratio": None, "p": None, "underpowered": True}
    return _fisher2x2(sub[col_a].map(pred_a).astype(bool), sub[col_b].map(pred_b).astype(bool))


def _concordance_perm(expected, actual, n_perm=20000, seed=0):
    """Agreement rate of expected-vs-actual class labels + a seeded label-permutation p."""
    import pandas as pd
    df = pd.DataFrame({"e": expected, "a": actual}).dropna()
    n = len(df)
    if n == 0:
        return {"n": 0, "agreement": None, "p_perm": None, "underpowered": True}
    e = df["e"].to_numpy(); a = df["a"].to_numpy()
    obs = float((e == a).mean())
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        if (rng.permutation(e) == a).mean() >= obs:
            ge += 1
    return {"n": n, "agreement": round(obs, 3), "p_perm": round((ge + 1) / (n_perm + 1), 5),
            "underpowered": n < UNDERPOWERED_N}


def validate(ctx) -> dict:
    t = build_table(ctx)
    anchored_mask = t["anchored_driver"].notna()

    # 1) anchored-driver ELN-expected risk vs the clinician's ELN_risk (3-class concordance)
    eln_conc = _concordance_perm(t["expected_eln"], t["eln_actual"])
    # 2) adverse anchored driver vs ELN Adverse (support = anchored & ELN)
    adv_vs_eln = _assoc(t, "expected_eln", lambda x: x == "Adverse", "eln_actual", lambda x: x == "Adverse")
    # 3) favorable anchored driver vs clinical responder (support = anchored & response)
    fav_vs_resp = _assoc(t, "expected_eln", lambda x: x == "Favorable", "responder", lambda x: x is True)
    # 4) LSC primitive-stem (p-LSC) vs ELN Adverse (support = confident-LSC & ELN) -- the most
    #    INDEPENDENT signal (RNA-derived stemness vs genetics-defined risk)
    lsc_vs_eln = _assoc(t, "lsc_class", lambda x: str(x).startswith("p"), "eln_actual", lambda x: x == "Adverse")

    # 5) survival by driver-risk (descriptive only; n is tiny)
    surv = t.dropna(subset=["survival"])
    surv_by_risk = {}
    for grp in ("Favorable", "Intermediate", "Adverse"):
        v = surv.loc[surv["expected_eln"] == grp, "survival"]
        if len(v):
            surv_by_risk[grp] = {"n": int(len(v)), "median_months": round(float(v.median()), 1)}

    report = {
        "mode": "retrospective_validation",
        "caveat": ("Retrospective ASSOCIATIONS on sparse, partially-overlapping clinical labels — "
                   "NOT trained predictors. n and an underpowered flag (<%d) accompany every test; "
                   "a null is a legitimate finding. ELN-expected is a single-driver proxy for a "
                   "multi-factor clinical schema. CIRCULARITY: ELN risk is itself defined from driver "
                   "genetics, so the anchored-driver ELN concordance + the adverse-driver->ELN test "
                   "SHARE logic with the label — they validate the engine's driver EXTRACTION + anchor "
                   "selection (a correctness/data-integrity check), NOT a novel prediction. The "
                   "p-LSC->ELN test is the more INDEPENDENT signal (an RNA-derived stemness phenotype "
                   "tracking genetics-defined risk). The clinical_response column is HETEROGENEOUSLY "
                   "encoded across cohorts (some use risk-style Favorable/Adverse labels, not response "
                   "tags), so the favorable-driver->responder result is SUSPECT (likely a labeling "
                   "artifact, not biology) and is not interpreted; that column needs curation."
                   % UNDERPOWERED_N),
        "test_independence": {
            "anchored_eln_concordance": "shares-logic (extraction/anchor correctness check)",
            "adverse_driver_vs_eln_adverse": "shares-logic (extraction/anchor correctness check)",
            "favorable_driver_vs_responder": "SUSPECT — heterogeneous clinical_response encoding; likely artifact, do not interpret",
            "pLSC_vs_eln_adverse": "independent (RNA-derived stemness vs genetic risk)",
            "survival_by_expected_risk": "independent but tiny-n (descriptive only)",
        },
        "coverage": {
            "n_samples": int(len(t)),
            "n_anchored_driver": int(anchored_mask.sum()),
            "n_eln": int(t["eln_actual"].notna().sum()),
            "n_responder": int(t["responder"].notna().sum()),
            "n_survival": int(t["survival"].notna().sum()),
            "n_lsc_confident": int(t["lsc_class"].notna().sum()),
        },
        "anchored_driver_distribution": {str(k): int(v) for k, v in
                                         t["anchored_driver"].dropna().value_counts().items()},
        "tests": {
            "anchored_eln_concordance": eln_conc,
            "adverse_driver_vs_eln_adverse": adv_vs_eln,
            "favorable_driver_vs_responder": fav_vs_resp,
            "pLSC_vs_eln_adverse": lsc_vs_eln,
            "survival_by_expected_risk": surv_by_risk,
        },
    }
    return report


def _write_md(ctx, report) -> str:
    c, cov, ts = report["caveat"], report["coverage"], report["tests"]
    L = ["# MOSAIC-AML retrospective clinical validation", "", f"_{c}_", "",
         "## Coverage (the binding constraint)",
         f"- samples: {cov['n_samples']}; anchored driver: {cov['n_anchored_driver']}; "
         f"ELN: {cov['n_eln']}; clinical response: {cov['n_responder']}; "
         f"survival: {cov['n_survival']}; confident LSC: {cov['n_lsc_confident']}", "",
         "## Tests"]

    indep = report.get("test_independence", {})

    def line(name, key, r, extra=""):
        tag = indep.get(key, "")
        tagstr = f" _[{tag}]_" if tag else ""
        if not r or r.get("n", 0) == 0:
            return f"- **{name}:** no overlapping data.{tagstr}"
        flag = " ⚠️ UNDERPOWERED" if r.get("underpowered") else ""
        if "agreement" in r:
            return (f"- **{name}:** agreement {r['agreement']} (n={r['n']}, "
                    f"perm p={r['p_perm']}){flag}{extra}{tagstr}")
        return (f"- **{name}:** OR={r.get('odds_ratio')} (n={r['n']}, Fisher p={r['p']}, "
                f"2x2={r['table']}){flag}{extra}{tagstr}")

    L.append(line("Anchored-driver ELN-expected vs clinician ELN (3-class concordance)",
                  "anchored_eln_concordance", ts["anchored_eln_concordance"]))
    L.append(line("Adverse anchored driver -> ELN Adverse",
                  "adverse_driver_vs_eln_adverse", ts["adverse_driver_vs_eln_adverse"]))
    L.append(line("Favorable anchored driver -> clinical responder",
                  "favorable_driver_vs_responder", ts["favorable_driver_vs_responder"]))
    L.append(line("Primitive-LSC (p-LSC) -> ELN Adverse",
                  "pLSC_vs_eln_adverse", ts["pLSC_vs_eln_adverse"]))
    sbr = ts["survival_by_expected_risk"]
    if sbr:
        L.append("- **Median survival (months) by expected risk** (descriptive, tiny n): "
                 + "; ".join(f"{k}: {v['median_months']} (n={v['n']})" for k, v in sbr.items()))
    L += ["", "### Anchored-driver distribution (cohort)",
          ", ".join(f"{k}={v}" for k, v in report["anchored_driver_distribution"].items())]
    fp = ctx.path("VALIDATION.md")
    with open(fp, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return fp


def run_validation(ctx) -> dict:
    report = validate(ctx)
    ctx.save_json(report, "validation_report.json")
    _write_md(ctx, report)
    return report
