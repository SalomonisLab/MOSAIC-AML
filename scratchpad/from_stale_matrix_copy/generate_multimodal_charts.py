#!/usr/bin/env python3
"""Generate PDF bar charts of Sensitivity, Specificity, and F1 for the MATRIX-AML
multimodal (all modalities) mutation predictor, from the sealed held-out reports.

Uses the EXISTING predict_* patient_report.json files in runs/ (which contain the
held-out predictions with true_label). These are the sealed, honest held-out scores
from the deployed late-fusion linSVM + 8-modality optimised weights system.

Produces:
  1) multimodal_all_mutations.pdf   — grouped bar chart (Sens/Spec/F1) with a dot
     per mutation overlaid, showing all mutations together.
  2) multimodal_per_mutation.pdf    — one page per mutation, each with 3 bars.
  3) combined_comparison.pdf        — side-by-side comparison: BEAT-AML RNA-Seq vs Multimodal.
"""

import os, sys, json, warnings, glob
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as mtick
from matplotlib.lines import Line2D

# ---- paths ----
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(os.path.dirname(HERE), "runs")
OUT_DIR = os.path.join(HERE, "charts")
os.makedirs(OUT_DIR, exist_ok=True)

# ---- load all held-out predict_* reports ----
print("Scanning held-out reports in:", RUNS)
all_preds = []   # list of {mutation, probability, call, threshold, true_label, sample}
for d in sorted(os.listdir(RUNS)):
    if not d.startswith("predict_"):
        continue
    rp = os.path.join(RUNS, d, "patient_report.json")
    if not os.path.exists(rp):
        continue
    with open(rp) as f:
        rep = json.load(f)
    sample = rep.get("sample_key", d)
    for p in rep.get("mutation_predictions", []):
        if p.get("true_label") is not None:
            all_preds.append({
                "sample": sample,
                "mutation": p["mutation"],
                "probability": p.get("probability"),
                "call": p.get("call"),
                "threshold": p.get("threshold", 0.5),
                "true_label": p["true_label"],
                "n_modalities": p.get("n_modalities", 0),
            })

df = pd.DataFrame(all_preds)
print(f"  {len(df)} predictions from {df['sample'].nunique()} held-out samples, "
      f"{df['mutation'].nunique()} mutations")

# ---- compute per-mutation metrics ----
results = {}
for mut, grp in df.groupby("mutation"):
    truth = (grp["true_label"] == "present").astype(int).values
    calls = (grp["call"] == "present").astype(int).values

    tp = int(((calls == 1) & (truth == 1)).sum())
    tn = int(((calls == 0) & (truth == 0)).sum())
    fp = int(((calls == 1) & (truth == 0)).sum())
    fn = int(((calls == 0) & (truth == 1)).sum())

    n_pos = int(truth.sum())
    n_neg = int((truth == 0).sum())

    # Skip mutations with 0 positives or 0 negatives (can't compute both sens & spec)
    if n_pos == 0 or n_neg == 0:
        print(f"  SKIP {mut}: n+={n_pos}, n-={n_neg} (no both classes)")
        continue

    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) > 0 else 0.0

    thr = grp["threshold"].iloc[0]

    results[mut] = {
        "sensitivity": round(sens, 3),
        "specificity": round(spec, 3),
        "f1": round(f1, 3),
        "precision": round(prec, 3),
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_total": n_pos + n_neg,
        "threshold": thr,
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
    }
    print(f"  {mut:25s}  sens={sens:.3f}  spec={spec:.3f}  F1={f1:.3f}  "
          f"(n+={n_pos}, n-={n_neg}, thr={thr})")

print(f"\n{len(results)} mutations with both classes in held-out")

# ---- save metrics ----
mdf = pd.DataFrame(results).T
mdf.index.name = "mutation"
mdf.to_csv(os.path.join(OUT_DIR, "multimodal_metrics.csv"))
print("Saved:", os.path.join(OUT_DIR, "multimodal_metrics.csv"))

# ---- styling ----
COL_SENS = "#2196F3"
COL_SPEC = "#4CAF50"
COL_F1   = "#FF9800"
DOT_COL  = "#263238"
BG_COL   = "#FAFAFA"
GRID_COL = "#E0E0E0"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.facecolor": BG_COL,
    "figure.facecolor": "white",
    "axes.edgecolor": "#BDBDBD",
    "axes.grid": True,
    "grid.color": GRID_COL,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.7,
})

# ============================================================
# PDF 1: All mutations — grouped bars + dots
# ============================================================
muts = sorted(results.keys(), key=lambda k: results[k]["f1"], reverse=True)
sens_vals = [results[m]["sensitivity"] for m in muts]
spec_vals = [results[m]["specificity"] for m in muts]
f1_vals   = [results[m]["f1"] for m in muts]

mean_sens = np.mean(sens_vals)
mean_spec = np.mean(spec_vals)
mean_f1   = np.mean(f1_vals)

fig, ax = plt.subplots(figsize=(10, 7))
x = np.arange(3)
bar_width = 0.55
bars = ax.bar(x, [mean_sens, mean_spec, mean_f1], bar_width,
              color=[COL_SENS, COL_SPEC, COL_F1], edgecolor="white", linewidth=1.2,
              alpha=0.85, zorder=3)

np.random.seed(42)
for i, vals in enumerate([sens_vals, spec_vals, f1_vals]):
    jitter = np.random.uniform(-0.15, 0.15, len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals,
               s=28, color=DOT_COL, alpha=0.55, zorder=5, edgecolors="white", linewidths=0.4)

for bar, val in zip(bars, [mean_sens, mean_spec, mean_f1]):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
            f"{val:.3f}", ha="center", va="bottom", fontweight="bold", fontsize=12, color="#333")

ax.set_xticks(x)
ax.set_xticklabels(["Sensitivity", "Specificity", "F1 Score"], fontsize=13, fontweight="bold")
ax.set_ylim(0, 1.12)
ax.set_ylabel("Score", fontsize=12)
ax.set_title("MATRIX-AML Multimodal (All 8 Modalities)\nAll Mutations — Sealed Held-Out — Mean ± Per-Mutation Dots",
             fontsize=14, fontweight="bold", pad=15)
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))

legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=DOT_COL, markersize=6,
           label=f"Individual mutations (n={len(muts)})", alpha=0.6),
    Line2D([0], [0], color=COL_SENS, lw=8, alpha=0.85, label="Mean Sensitivity"),
    Line2D([0], [0], color=COL_SPEC, lw=8, alpha=0.85, label="Mean Specificity"),
    Line2D([0], [0], color=COL_F1, lw=8, alpha=0.85, label="Mean F1"),
]
ax.legend(handles=legend_elements, loc="lower right", fontsize=9, framealpha=0.9)

ax.text(0.01, 0.99, f"Sealed held-out (n={df['sample'].nunique()} samples)\n"
        f"{len(muts)} mutations, 8 modalities\nlinSVM + optimised NNLS weights",
        transform=ax.transAxes, fontsize=8, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#ccc", alpha=0.9))

plt.tight_layout()
pdf1 = os.path.join(OUT_DIR, "multimodal_all_mutations.pdf")
with PdfPages(pdf1) as pp:
    pp.savefig(fig, dpi=300)
plt.close(fig)
print(f"Saved: {pdf1}")


# ============================================================
# PDF 2: Per-mutation bar charts
# ============================================================
pdf2 = os.path.join(OUT_DIR, "multimodal_per_mutation.pdf")
with PdfPages(pdf2) as pp:
    for mut in muts:
        r = results[mut]
        fig, ax = plt.subplots(figsize=(7, 5))
        vals = [r["sensitivity"], r["specificity"], r["f1"]]
        colors = [COL_SENS, COL_SPEC, COL_F1]
        labels = ["Sensitivity", "Specificity", "F1 Score"]
        xx = np.arange(3)

        bars = ax.bar(xx, vals, 0.55, color=colors, edgecolor="white", linewidth=1.2,
                       alpha=0.88, zorder=3)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                    f"{val:.3f}", ha="center", va="bottom", fontweight="bold", fontsize=13,
                    color="#333")

        ax.set_xticks(xx)
        ax.set_xticklabels(labels, fontsize=12, fontweight="bold")
        ax.set_ylim(0, 1.15)
        ax.set_ylabel("Score", fontsize=11)
        title_name = mut.replace("_", " ").replace("/", " / ")
        ax.set_title(f"MATRIX-AML Multimodal (All 8 Modalities)\n{title_name}",
                     fontsize=13, fontweight="bold", pad=12)
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))

        info = (f"n+ = {r['n_pos']}   n− = {r['n_neg']}   threshold = {r['threshold']}\n"
                f"TP={r['tp']}  FP={r['fp']}  FN={r['fn']}  TN={r['tn']}")
        ax.text(0.5, -0.18, info, transform=ax.transAxes, fontsize=8, ha="center",
                va="top", color="#555",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#ddd", alpha=0.9))

        plt.tight_layout(rect=[0, 0.08, 1, 1])
        pp.savefig(fig, dpi=300)
        plt.close(fig)

print(f"Saved: {pdf2}")


# ============================================================
# PDF 3: Side-by-side comparison (BEAT-AML RNA vs Multimodal)
# ============================================================
# Load the BEAT-AML RNA metrics for comparison
beataml_csv = os.path.join(OUT_DIR, "beataml_rna_metrics.csv")
if os.path.exists(beataml_csv):
    ba = pd.read_csv(beataml_csv, index_col=0)

    # Map BEAT-AML variant-level categories to gene-level for comparison
    # The multimodal system uses gene-level names, BEAT-AML uses variant-level
    # We'll map by taking the best variant-level performance per gene where applicable
    ba_gene_map = {
        "FLT3-ITD": ["FLT3_ITD"],
        "FLT3-TKD": ["FLT3_TKD_D835/I836"],
        "FLT3": ["FLT3_ITD", "FLT3_TKD_D835/I836", "FLT3_other_TKD_or_JM"],
        "NPM1": ["NPM1_exon12_frameshift"],
        "DNMT3A": ["DNMT3A_R882", "DNMT3A_nonR882"],
        "IDH1": ["IDH1_R132"],
        "IDH2": ["IDH2_R140", "IDH2_R172"],
        "NRAS": ["NRAS_G12", "NRAS_G13", "NRAS_Q61"],
        "TP53": ["TP53_hotspot_DBD", "TP53_LOF/splice/frameshift"],
        "RUNX1": ["RUNX1_LOF"],
        "WT1": ["WT1_LOF"],
        "TET2": ["TET2"],
        "ASXL1": ["ASXL1"],
        "SRSF2": ["SRSF2"],
        "STAG2": ["STAG2"],
        "BCOR": ["BCOR"],
        "CEBPA": ["CEBPA_bZIP", "CEBPA_Nterminal_frameshift/nonsense", "CEBPA_biallelic_or_double"],
        "KIT": ["KIT_D816"],
        "PTPN11": ["PTPN11_other"],
        "NF1": ["NF1"],
        "kmt2a": ["KMT2A_fusion"],
    }

    # For each multimodal mutation, find the best matching BEAT-AML entry
    comp = {}
    for mut in muts:
        mm = results[mut]
        # Try direct match, then gene-level lookup
        ba_variants = ba_gene_map.get(mut, [mut])
        ba_matches = [v for v in ba_variants if v in ba.index]
        if ba_matches:
            # Use the variant with the highest F1 as the representative
            best_v = max(ba_matches, key=lambda v: ba.loc[v, "f1"])
            comp[mut] = {
                "multimodal": mm,
                "beataml": {
                    "sensitivity": ba.loc[best_v, "sensitivity"],
                    "specificity": ba.loc[best_v, "specificity"],
                    "f1": ba.loc[best_v, "f1"],
                    "variant": best_v,
                }
            }

    if comp:
        comp_muts = sorted(comp.keys(), key=lambda k: comp[k]["multimodal"]["f1"], reverse=True)

        pdf3 = os.path.join(OUT_DIR, "combined_comparison.pdf")
        with PdfPages(pdf3) as pp:
            # Page 1: Summary comparison (grouped bars)
            fig, axes = plt.subplots(1, 3, figsize=(14, 6), sharey=True)
            metrics = ["sensitivity", "specificity", "f1"]
            metric_labels = ["Sensitivity", "Specificity", "F1 Score"]
            metric_colors = [COL_SENS, COL_SPEC, COL_F1]

            for ax, metric, mlabel, mcol in zip(axes, metrics, metric_labels, metric_colors):
                mm_vals = [comp[m]["multimodal"][metric] for m in comp_muts]
                ba_vals = [comp[m]["beataml"][metric] for m in comp_muts]
                mean_mm = np.mean(mm_vals)
                mean_ba = np.mean(ba_vals)

                x = np.arange(2)
                b = ax.bar(x, [mean_ba, mean_mm], 0.5,
                          color=["#78909C", mcol], edgecolor="white", linewidth=1.2,
                          alpha=0.85, zorder=3)

                np.random.seed(42)
                for i, vals in enumerate([ba_vals, mm_vals]):
                    jitter = np.random.uniform(-0.12, 0.12, len(vals))
                    ax.scatter(np.full(len(vals), i) + jitter, vals,
                               s=24, color=DOT_COL, alpha=0.45, zorder=5, edgecolors="white", linewidths=0.3)

                for bar, val in zip(b, [mean_ba, mean_mm]):
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                            f"{val:.3f}", ha="center", va="bottom", fontweight="bold", fontsize=11, color="#333")

                ax.set_xticks(x)
                ax.set_xticklabels(["BEAT-AML\nRNA-Seq", "Multimodal\n(8 mod.)"], fontsize=10, fontweight="bold")
                ax.set_ylim(0, 1.15)
                ax.set_title(mlabel, fontsize=13, fontweight="bold", color=mcol)
                ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))

            fig.suptitle("BEAT-AML RNA-Seq (No Imputation)  vs  MATRIX-AML Multimodal\n"
                         f"Matched Mutations (n={len(comp_muts)}) — Dots = Individual Mutations",
                         fontsize=14, fontweight="bold", y=1.02)

            legend_elements = [
                Line2D([0], [0], marker="o", color="w", markerfacecolor=DOT_COL, markersize=6,
                       label=f"Individual mutations (n={len(comp_muts)})", alpha=0.5),
            ]
            axes[1].legend(handles=legend_elements, loc="lower center", fontsize=8, framealpha=0.9)

            plt.tight_layout()
            pp.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)

            # Page 2: Per-mutation side-by-side comparison (horizontal)
            fig, ax = plt.subplots(figsize=(12, max(6, len(comp_muts) * 0.45)))

            y_pos = np.arange(len(comp_muts))
            bar_h = 0.35

            mm_f1 = [comp[m]["multimodal"]["f1"] for m in comp_muts]
            ba_f1 = [comp[m]["beataml"]["f1"] for m in comp_muts]

            ax.barh(y_pos + bar_h/2, mm_f1, bar_h, color=COL_F1, alpha=0.85,
                    label="Multimodal (8 mod.)", edgecolor="white", linewidth=0.8)
            ax.barh(y_pos - bar_h/2, ba_f1, bar_h, color="#78909C", alpha=0.85,
                    label="BEAT-AML RNA-Seq", edgecolor="white", linewidth=0.8)

            # Add value labels
            for i, (mm, ba) in enumerate(zip(mm_f1, ba_f1)):
                ax.text(mm + 0.01, i + bar_h/2, f"{mm:.2f}", va="center", fontsize=8, color="#333")
                ax.text(ba + 0.01, i - bar_h/2, f"{ba:.2f}", va="center", fontsize=8, color="#555")

            ax.set_yticks(y_pos)
            ax.set_yticklabels(comp_muts, fontsize=10)
            ax.set_xlabel("F1 Score", fontsize=12, fontweight="bold")
            ax.set_title("F1 Score Comparison by Mutation\nBEAT-AML RNA-Seq vs MATRIX-AML Multimodal",
                         fontsize=13, fontweight="bold")
            ax.set_xlim(0, 1.15)
            ax.legend(loc="lower right", fontsize=10, framealpha=0.9)
            ax.invert_yaxis()

            plt.tight_layout()
            pp.savefig(fig, dpi=300)
            plt.close(fig)

            # Pages 3+: Per-mutation comparison bar charts
            for mut in comp_muts:
                mm = comp[mut]["multimodal"]
                ba_r = comp[mut]["beataml"]

                fig, ax = plt.subplots(figsize=(8, 5))
                x = np.arange(3)
                w = 0.3

                bars1 = ax.bar(x - w/2, [ba_r["sensitivity"], ba_r["specificity"], ba_r["f1"]],
                              w, color="#78909C", alpha=0.85, label="BEAT-AML RNA-Seq",
                              edgecolor="white", linewidth=1)
                bars2 = ax.bar(x + w/2, [mm["sensitivity"], mm["specificity"], mm["f1"]],
                              w, color=[COL_SENS, COL_SPEC, COL_F1], alpha=0.85,
                              label="Multimodal (8 mod.)", edgecolor="white", linewidth=1)

                for bar, val in zip(bars1, [ba_r["sensitivity"], ba_r["specificity"], ba_r["f1"]]):
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                            f"{val:.3f}", ha="center", va="bottom", fontsize=10, color="#555")
                for bar, val in zip(bars2, [mm["sensitivity"], mm["specificity"], mm["f1"]]):
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                            f"{val:.3f}", ha="center", va="bottom", fontweight="bold", fontsize=10, color="#333")

                ax.set_xticks(x)
                ax.set_xticklabels(["Sensitivity", "Specificity", "F1 Score"],
                                   fontsize=12, fontweight="bold")
                ax.set_ylim(0, 1.18)
                ax.set_ylabel("Score", fontsize=11)
                title_name = mut.replace("_", " ").replace("/", " / ")
                ax.set_title(f"BEAT-AML RNA-Seq  vs  Multimodal\n{title_name}",
                             fontsize=13, fontweight="bold", pad=12)
                ax.legend(fontsize=9, loc="upper right", framealpha=0.9)
                ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))

                info = (f"Multimodal: n+={mm['n_pos']}  n−={mm['n_neg']}  TP={mm['tp']} FP={mm['fp']} FN={mm['fn']} TN={mm['tn']}\n"
                        f"BEAT-AML variant: {ba_r['variant']}")
                ax.text(0.5, -0.16, info, transform=ax.transAxes, fontsize=7.5, ha="center",
                        va="top", color="#555",
                        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#ddd", alpha=0.9))

                plt.tight_layout(rect=[0, 0.07, 1, 1])
                pp.savefig(fig, dpi=300)
                plt.close(fig)

        print(f"Saved: {pdf3}")
    else:
        print("No matching mutations between multimodal and BEAT-AML — skipping comparison PDF")
else:
    print("No BEAT-AML metrics found — skipping comparison PDF")

print(f"\nDONE. Charts in: {OUT_DIR}")
