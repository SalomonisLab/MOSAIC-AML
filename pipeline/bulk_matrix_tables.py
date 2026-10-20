#!/usr/bin/env python3
"""Write the bulk model's tables from deliverables/bulk_matrix.json.

  bulk_matrix.tsv          one row per (cohort, category)
  bulk_matrix_overall.tsv  one row per cohort: mean over categories and pooled over calls

bench_matrix.py calls write() after evaluating the model. Run on its own, this
rebuilds both tables from the committed JSON without the model or its data.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(os.path.dirname(HERE), "deliverables")

CATEGORY_FIELDS = ["n", "n_pos", "prevalence", "sensitivity", "specificity", "precision",
                   "f1", "auroc", "tp", "fp", "fn", "tn"]


def write(out, outd=OUTD):
    cohorts = out["cohorts"]
    with open(os.path.join(outd, "bulk_matrix.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("cohort\tcategory\tn\tn_pos\tprevalence\tsensitivity\tspecificity\tprecision\tf1\tauroc\tTP\tFP\tFN\tTN\n")
        for coh, c in cohorts.items():
            for cat, m in sorted(c["per_category"].items()):
                fh.write("\t".join(str(x) for x in [coh, cat] + [m[k] for k in CATEGORY_FIELDS]) + "\n")
    with open(os.path.join(outd, "bulk_matrix_overall.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("cohort\tn_categories\tmean_sensitivity\tmean_specificity\tmean_f1\tmean_auroc\t"
                 "pooled_sensitivity\tpooled_specificity\tpooled_TP\tpooled_FP\tpooled_FN\tpooled_TN\n")
        for coh, c in cohorts.items():
            s = c["overall"]
            if s:
                fh.write("%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\n" % (
                    coh, s["n_categories"], s["mean_sensitivity"], s["mean_specificity"],
                    s["mean_f1"], s["mean_auroc"], s["pooled_sensitivity"], s["pooled_specificity"],
                    s["pooled_TP"], s["pooled_FP"], s["pooled_FN"], s["pooled_TN"]))


if __name__ == "__main__":
    with open(os.path.join(OUTD, "bulk_matrix.json"), encoding="utf-8") as fh:
        write(json.load(fh))
    print("wrote deliverables/bulk_matrix.tsv + bulk_matrix_overall.tsv")
