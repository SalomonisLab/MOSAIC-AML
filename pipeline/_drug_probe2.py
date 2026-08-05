#!/usr/bin/env python3
"""Second survey pass: clinical covariates, wave/site metadata for leave-group-out validation, and how
the mutation table keys onto the drug-tested specimens."""
import os
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BA = os.path.join(ROOT, "data", "external", "beataml")

cl = pd.read_excel(os.path.join(BA, "beataml_wv1to4_clinical.xlsx"))
print("clinical", cl.shape)
print("columns:", list(cl.columns))
for c in cl.columns:
    if any(k in c.lower() for k in ("wave", "cohort", "site", "center", "institut", "dbgap", "specimen", "sample", "id")):
        v = cl[c]
        print("  %-34s nunique=%-6d ex=%s" % (c, v.nunique(), list(pd.Series(v.dropna().unique()).head(4))))

print("\n--- candidate outcome / covariate columns ---")
for c in cl.columns:
    if any(k in c.lower() for k in ("age", "sex", "gender", "vital", "overall", "survival", "eln", "risk",
                                    "blast", "wbc", "diagnos", "response", "treat", "priorMDS".lower(),
                                    "specimen_type", "type")):
        v = cl[c]
        print("  %-34s nunique=%-6d ex=%s" % (c, v.nunique(), list(pd.Series(v.dropna().astype(str).unique()).head(4))))

mp = pd.read_excel(os.path.join(BA, "beataml_waves1to4_sample_mapping.xlsx"))
print("\nsample mapping", mp.shape, list(mp.columns))
print(mp.head(4).to_string())

mut = pd.read_csv(os.path.join(BA, "mutations.txt"), sep="\t", low_memory=False, nrows=4000)
print("\nmutations cols:", list(mut.columns)[:30])
