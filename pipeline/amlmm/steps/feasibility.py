"""Step: feasibility — per-target trainability triage BEFORE any modeling.

For each candidate target: how many labeled samples, how many usable classes
(>= min_class_n), majority-class fraction, and cohort confounding (Cramer's V
between the label and the dataset of origin). Verdict yes/marginal/no. This is
the cheap analysis that stops you wasting effort on thinly labeled targets.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

from ..step import StepSpec, StepResult, register
from .. import targets as T


def _cramers_v(a, b):
    try:
        ct = pd.crosstab(a.astype(str), b.astype(str))
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            return None
        chi2 = chi2_contingency(ct)[0]
        n = ct.values.sum()
        k = min(ct.shape) - 1
        return float(np.sqrt(chi2 / (n * k))) if n and k else None
    except Exception:
        return None


def run(ctx, params) -> StepResult:
    s = ctx.tables["samples"]
    min_class_n = int(params["min_class_n"])
    rows = []
    for t, spec in T.TARGETS.items():
        lab = T.get_labels(ctx, t)
        labeled = lab.dropna()
        row = {"target": t, "kind": spec["kind"], "n_labeled": int(len(labeled))}
        if spec["kind"] == "numeric":
            row.update(n_classes_usable=None, classes_usable="", majority_frac=None,
                       cramers_v_vs_dataset=None,
                       verdict=("yes" if len(labeled) >= 60 else
                                "marginal" if len(labeled) >= 30 else "no"),
                       note="regression target (not classified by this pipeline)")
        else:
            usable = T.usable_classes(labeled, min_class_n)
            vc = labeled.value_counts()
            n_usable = int(vc.reindex(usable).sum()) if usable else 0
            maj = float(vc.iloc[0] / vc.sum()) if len(vc) else None
            cv = _cramers_v(labeled, s.loc[labeled.index, "dataset"])
            if len(usable) >= 2 and n_usable >= 40:
                verdict = "yes"
            elif len(usable) >= 2 and n_usable >= 20:
                verdict = "marginal"
            else:
                verdict = "no"
            note = ""
            if cv is not None and cv > 0.6:
                note = f"high cohort confounding (V={cv:.2f}) — use leave_one_cohort_out"
            row.update(n_classes_usable=len(usable),
                       classes_usable=";".join(map(str, usable)),
                       majority_frac=round(maj, 3) if maj is not None else None,
                       cramers_v_vs_dataset=round(cv, 3) if cv is not None else None,
                       verdict=verdict, note=note)
        rows.append(row)

    df = pd.DataFrame(rows).set_index("target")
    fp = ctx.save_table(df, "feasibility.tsv")
    ctx.set("feasibility", df.to_dict("index"))
    res = StepResult(name="feasibility", artifacts={"feasibility.tsv": fp})
    res.metrics = {
        "n_targets": int(len(df)),
        "trainable": [t for t, r in df.iterrows() if r["verdict"] == "yes"],
        "marginal": [t for t, r in df.iterrows() if r["verdict"] == "marginal"],
        "not_trainable": [t for t, r in df.iterrows() if r["verdict"] == "no"],
    }
    return res


STEP = register(StepSpec(
    name="feasibility",
    run=run,
    doc="Per-target trainability: labeled n, usable classes, cohort confounding, verdict.",
    params_schema={"min_class_n": {"default": 8, "doc": "min samples per class to count as usable"}},
    produces=["feasibility"],
))
