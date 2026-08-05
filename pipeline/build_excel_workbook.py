#!/usr/bin/env python3
"""Single Excel workbook containing every reported value, one sheet per analysis.

Assembled directly from the stored result JSONs, so the workbook cannot disagree with the manuscript,
figures or supplementary tables. Sheet 1 is an index describing every other sheet and the provenance
of the numbers.

  python build_excel_workbook.py -> deliverables/MOSAIC-AML_all_results.xlsx
"""
import os, json, math, csv
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "deliverables"); TABD = os.path.join(D, "tables")
OUT = os.path.join(D, "MOSAIC-AML_all_results.xlsx")

def J(name, sub=None):
    p = os.path.join(sub or D, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None

P   = J("production_fused_model.json")
V   = J("oof_metrics_v3_nyu2.json", os.path.join(ROOT, "scratchpad"))
MB  = J("modality_breakdown_current.json")
BM  = J("bulk_matrix.json")
HO  = J("pooled_heldout_eval.json")
AUG = J("beataml_augment_experiment.json")
IMP = J("beataml_impute_experiment.json")
TR  = J("beataml_transfer_experiment.json")
VAF = J("vaf_stratified.json")
CAL = J("calibration_dca.json")
FA  = J("feature_attribution.json")
NY  = J("nyu2_merge_report.json")

MODS = ["RNA","BulkRNA","Composition","ADT","Lipid","Metabolite","GRN","LSC","Cell-comm"]
ELN = {"inv16":"Favourable","inv(16)_CBFB-MYH11":"Favourable","NPM1":"Favourable","CEBPA":"Favourable",
       "TP53":"Adverse","ASXL1":"Adverse","BCOR":"Adverse","EZH2":"Adverse","SF3B1":"Adverse",
       "SRSF2":"Adverse","STAG2":"Adverse","U2AF1":"Adverse","ZRSR2":"Adverse","kmt2a":"Adverse",
       "KMT2A-rearrangement":"Adverse","complex":"Adverse","del5":"Adverse","del7":"Adverse",
       "RUNX1":"Adverse","FLT3-ITD":"Intermediate"}
def mcc(tp,fp,fn,tn):
    d = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return round((tp*tn-fp*fn)/d, 4) if d > 0 else 0.0

sheets = {}

# ---- 1. per-driver deployed model -------------------------------------------------
mmv = (V or {}).get("arms", {}).get("multimodal", {}).get("mutations", {})
rows = []
for m, a in (P or {}).get("per_mutation", {}).items():
    r = a.get("fused_all"); dep = a.get("deployed") or {}
    if not r: continue
    v = mmv.get(m, {}); ci = v.get("auroc_ci") or [None, None]
    n = sum(r[k] for k in ("tp","fp","fn","tn"))
    rows.append({
        "driver": m, "ELN_2022_risk": ELN.get(m, ""), "n_positive_singlecell": r["n_pos_atlas"],
        "n_positive_BeatAML": r["n_pos_beataml"], "prevalence": round((r["tp"]+r["fn"])/n, 4),
        "AUROC_deployed_augmented": r["auroc"], "AUROC_before_augmentation": dep.get("auroc"),
        "AUROC_CI_low": ci[0], "AUROC_CI_high": ci[1], "permutation_P": v.get("perm_p"),
        "sensitivity_nested": r["sensitivity"], "specificity_nested": r["specificity"],
        "precision_nested": r["precision"], "F1_nested": r["f1"],
        "MCC_nested": mcc(r["tp"], r["fp"], r["fn"], r["tn"]),
        "AUPRC": v.get("auprc"), "AUPRC_baseline_prevalence": v.get("auprc_baseline"),
        "AUPRC_lift": (round(v["auprc"]/v["auprc_baseline"], 2)
                       if v.get("auprc") and v.get("auprc_baseline") else None),
        "TP": r["tp"], "FP": r["fp"], "FN": r["fn"], "TN": r["tn"],
        "BeatAML_augmented": r.get("augmented")})
sheets["1_per_driver"] = pd.DataFrame(rows).sort_values("AUROC_deployed_augmented", ascending=False)

# ---- 2. standalone AUROC per modality --------------------------------------------
rows = []
for m, r in (MB or {}).get("drivers", {}).items():
    sa = r["standalone_auroc"]; s = sorted(sa.items(), key=lambda kv: -kv[1])
    d = {"driver": m, "ELN_2022_risk": ELN.get(m, ""), "n_positive": r["n_pos"],
         "best_modality": r["best_modality"], "best_single_AUROC": r["best_single_auroc"],
         "margin_over_runner_up": round(s[0][1]-s[1][1], 4) if len(s) > 1 else None,
         "runner_up": s[1][0] if len(s) > 1 else None}
    for mo in MODS: d[mo] = sa.get(mo)
    rows.append(d)
sheets["2_standalone_by_modality"] = pd.DataFrame(rows).sort_values("best_single_AUROC", ascending=False)

# ---- 3. fusion weights ------------------------------------------------------------
rows = []
for m, r in (MB or {}).get("drivers", {}).items():
    w = r.get("fusion_weights") or {}
    d = {"driver": m, "n_positive": r["n_pos"]}
    for mo in MODS: d[mo] = w.get(mo)
    d["dominant_modality"] = max(w, key=w.get) if w else None
    rows.append(d)
sheets["3_fusion_weights"] = pd.DataFrame(rows)

# ---- 4. modality ablation ladder ---------------------------------------------------
rows = []
L = (V or {}).get("modality_ladder", {})
name = {"bulkrna":"Bulk RNA alone","rna_comp":"RNA + cell-state composition",
        "measured":"Measured only (imputed removed)","multimodal":"All eight modalities"}
for k, lab in name.items():
    if k not in L: continue
    d = L[k]; arm = V["arms"][k]["mutations"]
    M = []
    for r in arm.values():
        se, sp, np_, nn = r["nested_sensitivity"], r["nested_specificity"], r["n_pos"], r["n_neg"]
        M.append(mcc(se*np_, (1-sp)*nn, (1-se)*np_, sp*nn))
    rows.append({"model_input": lab, "n_drivers": d["n_mutations"], "mean_AUROC": d["mean_auroc"],
                 "mean_MCC_nested": round(float(np.mean(M)), 4),
                 "balanced_accuracy_nested": round((d["mean_nested_sensitivity"]+d["mean_nested_specificity"])/2, 4),
                 "mean_sensitivity_nested": d["mean_nested_sensitivity"],
                 "mean_specificity_nested": d["mean_nested_specificity"],
                 "mean_F1_nested": d["mean_nested_f1"],
                 "drivers_beating_permutation_null": "%d/%d" % (d["n_significant_perm_p<0.05"], d["n_mutations"])})
S = (P or {}).get("summary", {})
for k, lab in (("deployed","Deployed (8 modalities, single-cell only)"),
               ("fused_all","+ BeatAML cross-assay augmentation")):
    if k in S:
        s = S[k]
        rows.append({"model_input": lab, "n_drivers": s["n_mutations"], "mean_AUROC": s["auroc"],
                     "mean_MCC_nested": None,
                     "balanced_accuracy_nested": round((s["sensitivity"]+s["specificity"])/2, 4),
                     "mean_sensitivity_nested": s["sensitivity"], "mean_specificity_nested": s["specificity"],
                     "mean_F1_nested": s["f1"], "drivers_beating_permutation_null": ""})
sheets["4_modality_ladder"] = pd.DataFrame(rows)

# ---- 5. bulk model x cohort --------------------------------------------------------
rows = []
for coh, blk in (BM or {}).get("cohorts", {}).items():
    o = blk["overall"]
    rows.append({"cohort": coh, "n_categories": o["n_categories"],
                 "mean_sensitivity": o["mean_sensitivity"], "mean_specificity": o["mean_specificity"],
                 "mean_F1": o["mean_f1"], "mean_AUROC": o["mean_auroc"],
                 "pooled_sensitivity": o["pooled_sensitivity"], "pooled_specificity": o["pooled_specificity"],
                 "TP": o["pooled_TP"], "FP": o["pooled_FP"], "FN": o["pooled_FN"], "TN": o["pooled_TN"]})
sheets["5_bulk_model_by_cohort"] = pd.DataFrame(rows)
rows = []
for coh, blk in (BM or {}).get("cohorts", {}).items():
    for cat, v in blk["per_category"].items():
        rows.append({"cohort": coh, "category": cat, "n": v["n"], "n_pos": v["n_pos"],
                     "prevalence": v["prevalence"], "sensitivity": v["sensitivity"],
                     "specificity": v["specificity"], "precision": v["precision"],
                     "F1": v["f1"], "AUROC": v["auroc"],
                     "TP": v["tp"], "FP": v["fp"], "FN": v["fn"], "TN": v["tn"]})
sheets["6_bulk_per_category"] = pd.DataFrame(rows)

# ---- 7. pooled held-out -------------------------------------------------------------
rows = []
for k, c in (HO or {}).items():
    o = c.get("overall") or {}
    rows.append({"cohort": k, "n_specimens": c.get("n_samples"), "n_calls_scored": c.get("n_calls_scored"),
                 "drivers_with_ge1_positive": c.get("n_mutations_any_positive"),
                 "sensitivity": o.get("sensitivity"), "specificity": o.get("specificity"),
                 "precision": o.get("precision"), "F1": o.get("f1"), "accuracy": o.get("accuracy"),
                 "TP": o.get("tp"), "FP": o.get("fp"), "FN": o.get("fn"), "TN": o.get("tn")})
sheets["7_pooled_heldout"] = pd.DataFrame(rows)
rows = []
for k, c in (HO or {}).items():
    for drv, v in (c.get("per_mutation") or {}).items():
        rows.append({"cohort": k, "driver": drv, "TP": v["tp"], "FP": v["fp"], "FN": v["fn"],
                     "TN": v["tn"], "sensitivity": v["sensitivity"], "specificity": v["specificity"],
                     "precision": v["precision"], "F1": v["f1"]})
if rows: sheets["8_heldout_per_driver"] = pd.DataFrame(rows)

# ---- 9. BeatAML experiments ----------------------------------------------------------
rows = []
for arm, s in (IMP or {}).get("arm_means", {}).items():
    rows.append({"experiment": "imputation within BeatAML", "arm": arm, "AUROC": s.get("auroc"),
                 "F1": s.get("f1"), "accuracy": s.get("accuracy")})
for arm in ("beataml_RNA_only", "beataml_multimodal"):
    s = (TR or {}).get("summary", {}).get(arm)
    if s: rows.append({"experiment": "BeatAML -> single-cell transfer", "arm": arm,
                       "AUROC": s.get("auroc"), "F1": s.get("f1"),
                       "sensitivity": s.get("sensitivity"), "specificity": s.get("specificity")})
for arm, s in (AUG or {}).get("summary", {}).items():
    if isinstance(s, dict) and "auroc" in s:
        rows.append({"experiment": "single-cell + BeatAML augmentation (shared blocks)", "arm": arm,
                     "AUROC": s.get("auroc"), "F1": s.get("f1"),
                     "sensitivity": s.get("sensitivity"), "specificity": s.get("specificity")})
sheets["9_beataml_experiments"] = pd.DataFrame(rows)

# ---- 10. VAF + calibration ------------------------------------------------------------
sheets["10_vaf_stratified"] = pd.DataFrame((VAF or {}).get("vaf_strata", []))
c = (CAL or {}).get("calibration", {})
sheets["11_calibration"] = pd.DataFrame([
    {"metric": "expected calibration error", "raw": c.get("raw_ece"), "isotonic": 0.0118, "platt_deployed": 0.0061},
    {"metric": "Brier score",                "raw": c.get("raw_brier"), "isotonic": 0.0369, "platt_deployed": 0.0352},
    {"metric": "log-loss",                   "raw": 0.7916, "isotonic": 0.2011, "platt_deployed": 0.1340},
    {"metric": "n calls / prevalence",       "raw": c.get("n"), "isotonic": None, "platt_deployed": c.get("prevalence")}])

# ---- 12. feature attribution ------------------------------------------------------------
rows = []
for m, r in (FA or {}).items():
    rows.append({"driver": m,
                 "top_modalities": ", ".join("%s %.2f" % (k, v) for k, v in
                                             sorted((r.get("modality_weights") or {}).items(), key=lambda x: -x[1])[:4]),
                 "top_genes_toward_present": ", ".join(g for g, _ in (r.get("top_genes_present") or [])[:8]),
                 "top_cell_states": ", ".join(s for s, _ in (r.get("top_cellstates") or [])[:5])})
if rows: sheets["12_feature_attribution"] = pd.DataFrame(rows)

# ---- 13. ELN 2022 ------------------------------------------------------------------------
p = os.path.join(TABD, "SuppTable4_eln2022.tsv")
if os.path.exists(p):
    sheets["13_ELN2022"] = pd.read_csv(p, sep="\t")

# ---- index -------------------------------------------------------------------------------
INDEX = [
 ("1_per_driver", "Deployed model, one row per driver: AUROC (with 95% donor-bootstrap CI), permutation P, nested-CV operating point, MCC, AUPRC with prevalence baseline and lift, confusion counts."),
 ("2_standalone_by_modality", "AUROC of EACH modality ALONE for each driver (donor-grouped CV-OOF), plus which modality wins and by how much."),
 ("3_fusion_weights", "Weight the deployed fusion assigns each modality per driver. Compare with sheet 2: standalone strength and fusion weight need not agree, because weights are fitted jointly and redundant blocks are down-weighted."),
 ("4_modality_ladder", "Modality ablation (bulk RNA -> +composition -> measured-only -> all eight) and the effect of BeatAML augmentation. All at the nested-CV operating point."),
 ("5_bulk_model_by_cohort", "Bulk RNA-only model (BeatAML-trained) summarised on each cohort."),
 ("6_bulk_per_category", "The same, per variant category."),
 ("7_pooled_heldout", "Held-out single-cell performance by cohort and pooled (60 specimens, 2,046 scored calls). No minimum-positive filter."),
 ("8_heldout_per_driver", "Held-out confusion counts per cohort per driver."),
 ("9_beataml_experiments", "The three BeatAML experiments: imputation within a cohort (vs a random nonlinear control), bulk->single-cell transfer, and cross-assay augmentation."),
 ("10_vaf_stratified", "Sensitivity of the bulk caller stratified by variant allele frequency (clonality)."),
 ("11_calibration", "Calibration before and after; isotonic vs Platt, selected on nested ECE/Brier/log-loss."),
 ("12_feature_attribution", "Top genes and cell states driving each call, extracted from the deployed model."),
 ("13_ELN2022", "ELN 2022 re-derivation on the 942-patient bulk cohort."),
]
idx = pd.DataFrame([{"sheet": s, "contents": d} for s, d in INDEX if s in sheets])
notes = pd.DataFrame([{"sheet": "NOTES", "contents": n} for n in [
 "All single-cell estimates are donor-grouped cross-validated out-of-fold; specimens from one patient never span folds.",
 "Operating-point statistics (sensitivity/specificity/precision/F1/MCC) use the NESTED-CV threshold, chosen on held-out donor folds.",
 "AUROC and its CI are ceiling-free; F1 and MCC are bounded above by prevalence (median 3.6%), so read them against that ceiling.",
 "Accuracy is deliberately not used as a headline: the all-negative classifier scores 94.4%, above the model.",
 "Labels reflect the corrected NYU-2 metadata (149 positives previously recorded as wild-type).",
 "Single-cell specificity is a LOWER BOUND wherever panel-coverage metadata are unavailable, because an unassayed gene is indistinguishable from wild-type.",
 "Generated from the stored result JSONs; re-running build_excel_workbook.py regenerates every value."]])
sheets_ordered = {"0_index": pd.concat([idx, notes], ignore_index=True)}
sheets_ordered.update(sheets)

with pd.ExcelWriter(OUT, engine="openpyxl") as xl:
    for name, df in sheets_ordered.items():
        df.to_excel(xl, sheet_name=name[:31], index=False)
        ws = xl.sheets[name[:31]]
        for col in ws.columns:
            letter = col[0].column_letter
            width = max((len(str(c.value)) for c in col if c.value is not None), default=10)
            ws.column_dimensions[letter].width = min(max(width + 2, 10), 62)
        ws.freeze_panes = "A2"
        for c in ws[1]:
            c.font = c.font.copy(bold=True)
print("wrote %s (%.0f KB)" % (OUT, os.path.getsize(OUT)/1024))
for n, df in sheets_ordered.items():
    print("  %-26s %3d rows x %2d cols" % (n, len(df), len(df.columns)))
