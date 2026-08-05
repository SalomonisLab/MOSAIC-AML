#!/usr/bin/env python3
"""Survey the BeatAML2 ex-vivo inhibitor panel before any modelling: what is measurable, on how many
patients, with what curve quality, and how much of it overlaps the expression matrix we already hold."""
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BA = os.path.join(ROOT, "data", "external", "beataml")

d = pd.read_csv(os.path.join(BA, "beataml_probit_curve_fits_v4_dbgap.txt"), sep="\t", low_memory=False)
print("probit rows %d  cols %d" % d.shape)
print("type:", d["type"].value_counts().to_dict())
print("status:", d["status"].value_counts().to_dict())
print("paper_inclusion:", d["paper_inclusion"].value_counts().to_dict())
print("converged:", d["converged"].value_counts().to_dict())
print("curve_type:", d["curve_type"].value_counts().to_dict())
print("all_gt_50:", d["all_gt_50"].value_counts().to_dict(), " all_lt_50:", d["all_lt_50"].value_counts().to_dict())
print("subjects %d | rnaseq samples %d (non-null %d) | inhibitors %d"
      % (d["dbgap_subject_id"].nunique(), d["dbgap_rnaseq_sample"].nunique(),
         d["dbgap_rnaseq_sample"].notna().sum(), d["inhibitor"].nunique()))
print("auc range %.1f - %.1f  median %.1f" % (d["auc"].min(), d["auc"].max(), d["auc"].median()))
print("auc NaN:", int(d["auc"].isna().sum()))

sa = d[d["type"] == "single-agent"]
print("\nsingle-agent rows %d | inhibitors %d" % (len(sa), sa["inhibitor"].nunique()))
cnt = sa.groupby("inhibitor")["dbgap_subject_id"].nunique().sort_values(ascending=False)
print("patients-per-drug: median %d  max %d  min %d" % (cnt.median(), cnt.max(), cnt.min()))
for k in (30, 50, 100, 200, 300):
    print("  drugs with >= %3d patients: %d" % (k, int((cnt >= k).sum())))
print("\ntop 15 drugs by n:", [(i, int(v)) for i, v in cnt.head(15).items()])

# ---- expression overlap ----
hdr = pd.read_csv(os.path.join(BA, "norm_exp.txt"), sep="\t", nrows=1)
exp_cols = set(c for c in hdr.columns if c.startswith("BA"))
print("\nnorm_exp sample columns: %d  (e.g. %s)" % (len(exp_cols), list(exp_cols)[:3]))
have = sa["dbgap_rnaseq_sample"].dropna().astype(str)
print("single-agent rows with an rnaseq sample present in norm_exp: %d / %d" %
      (int(have.isin(exp_cols).sum()), len(sa)))
ok = sa[sa["dbgap_rnaseq_sample"].astype(str).isin(exp_cols)]
print("=> usable (drug x expressed-specimen) pairs: %d over %d specimens, %d subjects, %d drugs"
      % (len(ok), ok["dbgap_rnaseq_sample"].nunique(), ok["dbgap_subject_id"].nunique(),
         ok["inhibitor"].nunique()))
c2 = ok.groupby("inhibitor")["dbgap_rnaseq_sample"].nunique().sort_values(ascending=False)
for k in (30, 50, 100, 150, 200, 300):
    print("   with expression, drugs with >= %3d specimens: %d" % (k, int((c2 >= k).sum())))

# serial specimens per subject (leakage risk)
per = ok.groupby("dbgap_subject_id")["dbgap_rnaseq_sample"].nunique()
print("\nsubjects with >1 expressed specimen: %d / %d (max %d)"
      % (int((per > 1).sum()), len(per), int(per.max())))

fam = pd.read_excel(os.path.join(BA, "beataml_drug_families.xlsx"))
print("\ndrug families file: %s  cols=%s" % (fam.shape, list(fam.columns)[:10]))
print(fam.head(6).to_string())
