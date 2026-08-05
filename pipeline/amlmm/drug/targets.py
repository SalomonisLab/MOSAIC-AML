"""Inhibitor annotation: what each compound hits, how it is used clinically, and how it fails.

A curated table (`knowledge/drug_annotation.tsv`) rather than code, so the lab can extend it without
touching the pipeline -- same pattern as `biomarker_drug.tsv`. Every column exists to stop a specific
failure mode of "rank the drugs by predicted sensitivity":

  targets / family / mechanism  the mechanistic model needs a gene set to interrogate
  clinical_tier                 an exquisitely sensitive tool compound must never outrank an approved
                                agent in the same list -- rankings are produced PER TIER
  analogue                      a research compound whose clinical stand-in exists is far more useful
                                than one without (AGI-6780 -> enasidenib); the report says so
  exposure                      whether the assay's concentration window plausibly overlaps human
                                exposure -- coarse three-level curation, NOT a PK model
  resistance                    known bypass routes, checked against the patient's own features

`exposure` is deliberately coarse. The quantitative half of that question is answered from data, not
curation: `assay_reliance()` asks whether a specimen's sensitivity only appears at the top of the
tested concentration range, which is the signature of a non-specific effect.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
KB = os.path.join(os.path.dirname(HERE), "knowledge", "drug_annotation.tsv")

TIERS = ["approved_AML", "approved_other", "trial", "research"]
TIER_LABEL = {"approved_AML": "approved in AML",
              "approved_other": "approved, other indication (off-label in AML)",
              "trial": "clinical-trial agent",
              "research": "research-only compound"}
# clinical-evidence term of the utility score; an approved AML drug is worth more than an equally
# sensitive tool compound because it can actually be given to the patient
TIER_SCORE = {"approved_AML": 1.00, "approved_other": 0.70, "trial": 0.45, "research": 0.15}
EXPOSURE_SCORE = {"clinically_achievable": 1.0, "borderline": 0.6, "not_established": 0.25}

# The curated `family` is fine-grained (63 families over 118 drugs, many singletons) -- too fine for a
# hierarchical model to borrow strength across. `family_group` is the coarser pathway level at which
# pooling is biologically defensible: two VEGFR inhibitors really should share a prior, an HDAC
# inhibitor and a MEK inhibitor should not.
FAMILY_GROUP = {
    "VEGFR": "RTK", "VEGFR/FGFR": "RTK", "VEGFR/KIT": "RTK", "VEGFR/RET": "RTK", "VEGFR/RAF": "RTK",
    "VEGFR/FLT3": "RTK", "FGFR": "RTK", "MET": "RTK", "MET/VEGFR": "RTK", "KIT": "RTK",
    "ABL/KIT": "RTK", "CSF1R": "RTK", "IGF1R": "RTK", "ALK": "RTK", "ALK/MET": "RTK",
    "NTRK/ALK": "RTK", "ERBB": "RTK", "EGFR": "RTK", "ABL/SRC": "RTK", "SRC": "RTK",
    "FLT3": "FLT3", "FLT3/RAF": "FLT3", "FLT3/Aurora": "FLT3", "ABL/FLT3": "FLT3",
    "MEK": "MAPK", "RAF": "MAPK",
    "PI3K": "PI3K_AKT_mTOR", "PI3K/mTOR": "PI3K_AKT_mTOR", "AKT": "PI3K_AKT_mTOR",
    "AKT/PDK1": "PI3K_AKT_mTOR", "mTOR": "PI3K_AKT_mTOR",
    "JAK": "JAK_STAT", "STAT": "JAK_STAT",
    "CDK": "cell_cycle", "Aurora": "cell_cycle", "Aurora/CDK": "cell_cycle", "PLK": "cell_cycle",
    "DDR": "cell_cycle",
    "BET": "epigenetic", "HDAC": "epigenetic", "hypomethylating": "epigenetic", "IMiD": "epigenetic",
    "BCL2": "apoptosis", "IAP": "apoptosis", "p53": "apoptosis",
    "proteasome": "proteostasis", "HSP90": "proteostasis", "nuclear export": "proteostasis",
    "BTK": "immune_signalling", "SYK": "immune_signalling", "PKC": "immune_signalling",
    "NF-kB": "NF-kB", "p38": "stress_MAPK", "PKA": "stress_MAPK", "CAMKK": "stress_MAPK",
    "IDH": "metabolic", "metabolic": "metabolic", "nuclear receptor": "metabolic",
    "antimetabolite": "chemotherapy",
    "NOTCH": "developmental", "WNT": "developmental", "TGFb": "developmental",
}

_CACHE = {}


def annotation(path=KB):
    """inhibitor -> annotation dict (targets as a list; empty fields normalised to None)."""
    if path in _CACHE:
        return _CACHE[path]
    d = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    out = {}
    for _, r in d.iterrows():
        def _split(v):
            v = (v or "").strip()
            return [] if v in ("", "-", "—") else [x.strip() for x in v.split("|") if x.strip()]
        out[r["inhibitor"].strip()] = {
            "inhibitor": r["inhibitor"].strip(),
            "targets": _split(r["targets"]),
            "family": (r["family"] or "").strip() or "other",
            "mechanism": (r["mechanism"] or "").strip(),
            "clinical_tier": (r["clinical_tier"] or "research").strip(),
            "analogue": _split(r.get("analogue", "")),
            "exposure": (r["exposure"] or "not_established").strip(),
            "resistance": _split(r.get("resistance", "")),
        }
        out[r["inhibitor"].strip()]["family_group"] = FAMILY_GROUP.get(
            out[r["inhibitor"].strip()]["family"], "other")
    _CACHE[path] = out
    return out


def get(inhibitor):
    """Annotation for one inhibitor, or a permissive 'unknown' record so nothing crashes on a new drug."""
    a = annotation().get(str(inhibitor).strip())
    if a is not None:
        return a
    return {"inhibitor": inhibitor, "targets": [], "family": "other", "family_group": "other",
            "mechanism": "", "clinical_tier": "research", "analogue": [],
            "exposure": "not_established", "resistance": [], "unannotated": True}


def families(inhibitors=None, level="family_group"):
    """{group -> [inhibitors]}. `level='family_group'` is the pooling level for the hierarchical model."""
    ann = annotation()
    keys = list(ann) if inhibitors is None else [i for i in inhibitors if i in ann]
    out = {}
    for k in keys:
        out.setdefault(ann[k][level], []).append(k)
    return out


def clinical_score(inhibitor):
    a = get(inhibitor)
    s = TIER_SCORE.get(a["clinical_tier"], 0.15)
    if a["clinical_tier"] == "research" and a["analogue"]:
        s = max(s, 0.35)                       # a research compound with a clinical stand-in is actionable
    return s


def exposure_score(inhibitor):
    return EXPOSURE_SCORE.get(get(inhibitor)["exposure"], 0.25)


def assay_reliance(curves, inhibitor):
    """Does this drug's *measured* sensitivity live at the top of the tested concentration range?

    Computed from the curve fits, not curated. For each drug we take the fraction of sensitive-tail
    specimens whose fitted IC50 sits in the top decade of the tested window. A high fraction means the
    apparent sensitivity is only reachable at concentrations where selectivity is doubtful, and the
    pharmacology agent should discount the recommendation regardless of how confident the model is.
    """
    g = curves[curves["inhibitor"] == inhibitor]
    if not len(g):
        return None
    ic50 = pd.to_numeric(g["ic50"], errors="coerce")
    mx = pd.to_numeric(g["max_conc"], errors="coerce")
    sens = g["resp_class"].eq("sensitive")
    hi = (ic50 >= mx / 10.0)
    frac_all = float(np.nanmean(hi.astype(float))) if len(g) else np.nan
    frac_sens = float(np.nanmean(hi[sens].astype(float))) if int(sens.sum()) else np.nan
    return {"frac_ic50_top_decade": round(frac_all, 3),
            "frac_sensitive_ic50_top_decade": None if frac_sens != frac_sens else round(frac_sens, 3),
            "max_conc_uM": None if mx.dropna().empty else float(mx.dropna().mode().iloc[0]),
            "min_conc_uM": None if pd.to_numeric(g["min_conc"], errors="coerce").dropna().empty
                           else float(pd.to_numeric(g["min_conc"], errors="coerce").dropna().mode().iloc[0])}


if __name__ == "__main__":
    ann = annotation()
    print("annotated inhibitors: %d" % len(ann))
    from collections import Counter
    print("clinical tiers:", Counter(a["clinical_tier"] for a in ann.values()))
    print("exposure:", Counter(a["exposure"] for a in ann.values()))
    fam = families()
    print("families: %d" % len(fam))
    for f, v in sorted(fam.items(), key=lambda kv: -len(kv[1]))[:12]:
        print("  %-16s %2d  %s" % (f, len(v), ", ".join(sorted(v)[:5]) + ("..." if len(v) > 5 else "")))
    n_no_target = [k for k, a in ann.items() if not a["targets"]]
    print("without a gene-level target:", n_no_target)
