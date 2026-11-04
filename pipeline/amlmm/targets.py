"""Prediction targets + label extraction (shared by feasibility and assembly).

Defines the candidate targets and a single source of truth for turning the
canonical tables into a sample-level label vector, so feasibility scoring and
feature assembly never disagree about what a label means.
"""
from __future__ import annotations
import pandas as pd
from .hooks import NON_DRIVER_LABELS

TARGETS = {
    "subtype":            {"kind": "multiclass", "source": "annotation",
                           "desc": "genetic driver (canonicalized; non-driver labels dropped)"},
    "disease_category":   {"kind": "multiclass", "source": "disease_category"},
    "WHO_classification": {"kind": "multiclass", "source": "WHO_classification"},
    "FAB":                {"kind": "multiclass", "source": "FAB"},
    "ELN_risk":           {"kind": "multiclass", "source": "ELN_risk"},
    "sex":                {"kind": "binary",     "source": "sex"},
    "is_pediatric":       {"kind": "binary",     "source": "is_pediatric"},
    "age":                {"kind": "numeric",    "source": "age"},
}


def get_labels(ctx, target) -> pd.Series:
    """Sample-level label series (index = sample_key), missing -> NaN."""
    s = ctx.tables["samples"]
    if target == "subtype":
        m = ctx.hooks.canonical_label_map(s["annotation"].tolist())
        lab = s["annotation"].map(m)
        return lab.where(~lab.isin(NON_DRIVER_LABELS))
    spec = TARGETS.get(target)
    if spec is None or spec["source"] not in s.columns:
        return pd.Series(index=s.index, dtype=object)
    lab = s[spec["source"]]
    bad = lab.isna() | (lab.astype(str).str.lower() == "nan") | (lab.astype(str) == "")
    return lab.where(~bad)


def usable_classes(labels, min_class_n) -> list:
    vc = labels.dropna().value_counts()
    return list(vc[vc >= int(min_class_n)].index)
