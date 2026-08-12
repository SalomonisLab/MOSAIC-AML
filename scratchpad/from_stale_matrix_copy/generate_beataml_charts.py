#!/usr/bin/env python3
"""Generate PDF bar charts of Sensitivity, Specificity, and F1 for the BEAT-AML RNA-Seq
(no imputation) bulk mutation predictor.

Produces:
  1) beataml_rna_all_mutations.pdf   — one grouped bar chart (Sens/Spec/F1) with a dot
     per mutation overlaid on each bar, showing all mutations together.
  2) beataml_rna_per_mutation.pdf    — one page per mutation, each with a 3-bar chart
     (Sensitivity, Specificity, F1).

Method:  Replicates the exact 5-fold stratified CV-OOF from train_bulk_predictor.py
         (logL2, top-2500 variance, BeatAML training data), then applies the F1-max
         threshold stored in bulk_model_card.json to binarize calls and compute
         Sensitivity (recall), Specificity, and F1.
"""

import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as mtick

# ---- paths ----
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BUNDLE = r"C:\Users\krog5w\.gemini\antigravity\scratch\aml-bakeoff\bundle_data.npz"
CARD = os.path.join(HERE, "bulk_model_card.json")
OUT_DIR = os.path.join(HERE, "charts")
os.makedirs(OUT_DIR, exist_ok=True)

# ---- load ----
print("Loading bundle:", BUNDLE)
d = np.load(BUNDLE, allow_pickle=True)
genes = [str(g) for g in d["genes"]]
cats = [str(c) for c in d["drivers"]]
baX = d["ba_X"].astype(float)
baL = d["ba_L"].astype(float)

with open(CARD) as f:
    card = json.load(f)

# ---- replicate CV-OOF (exact same recipe as train_bulk_predictor.py) ----
VARCAP = 2500
MIN_POS = 6

def clog_z(X):
    cl = np.log2(np.clip(X, 0, None) + 1.0)
    mu = cl.mean(0); sd = cl.std(0); sd[sd == 0] = 1.0
    return (cl - mu) / sd

Zba = clog_z(baX)

def topvar(Z, cap=VARCAP):
    return np.argsort(Z.var(0))[::-1][:cap]

def pct(sorted_train, v):
    n = len(sorted_train); v = np.asarray(v, float); out = np.full(len(v), np.nan)
    ok = ~np.isnan(v)
    if n >= 2:
        out[ok] = np.searchsorted(sorted_train, v[ok], side="right") / n
    elif n == 1:
        out[ok] = 0.5
    return out

results = {}   # cat -> {sensitivity, specificity, f1, n_pos, n_neg, threshold, cv_auroc}

for j, cat in enumerate(cats):
    y = baL[:, j]
    ok = ~np.isnan(y)
    yv = y[ok].astype(int)
    Z = Zba[ok]

    if yv.sum() < MIN_POS or (yv == 0).sum() < MIN_POS:
        continue
    if cat not in card.get("categories", {}):
        continue

    # 5-fold CV-OOF decision scores
    oof = np.full(len(yv), np.nan)
    ns = min(5, int(yv.sum()), int((yv == 0).sum()))
    for tri, tei in StratifiedKFold(ns, shuffle=True, random_state=0).split(Z, yv):
        sel = topvar(Z[tri])
        est = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
        est.fit(Z[tri][:, sel], yv[tri])
        oof[tei] = est.decision_function(Z[tei][:, sel])

    om = ~np.isnan(oof)
    if om.sum() < 10:
        continue

    # orientation
    from sklearn.metrics import roc_auc_score
    try:
        a = roc_auc_score(yv[om], oof[om])
    except Exception:
        a = np.nan
    sign = 1.0 if (a != a or a >= 0.5) else -1.0
    cv_auroc = round(float(max(a, 1 - a)), 3) if a == a else None

    # percentile-calibrated probabilities
    sorted_oof = np.sort(oof[om])
    p = pct(sorted_oof, oof[om])
    if sign < 0:
        p = 1.0 - p

    # get stored threshold from model card
    thr = card["categories"][cat].get("threshold", 0.5)

    # binary calls
    calls = (p >= thr).astype(int)
    truth = yv[om].astype(int)

    tp = int(((calls == 1) & (truth == 1)).sum())
    tn = int(((calls == 0) & (truth == 0)).sum())
    fp = int(((calls == 1) & (truth == 0)).sum())
    fn = int(((calls == 0) & (truth == 1)).sum())

    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) > 0 else 0.0

    results[cat] = {
        "sensitivity": round(sens, 3),
        "specificity": round(spec, 3),
        "f1": round(f1, 3),
        "precision": round(prec, 3),
        "cv_auroc": cv_auroc,
        "n_pos": int(yv.sum()),
        "n_neg": int((yv == 0).sum()),
        "threshold": thr,
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
    }
    print(f"  {cat:40s}  sens={sens:.3f}  spec={spec:.3f}  F1={f1:.3f}  "
          f"(n+={int(yv.sum())}, thr={thr:.3f}, AUC={cv_auroc})")

print(f"\n{len(results)} mutations computed")

# ---- save metrics table ----
df = pd.DataFrame(results).T
df.index.name = "mutation"
df.to_csv(os.path.join(OUT_DIR, "beataml_rna_metrics.csv"))
print("Saved metrics CSV to", os.path.join(OUT_DIR, "beataml_rna_metrics.csv"))

# ---- styling ----
# Professional color palette
COL_SENS = "#2196F3"   # blue
COL_SPEC = "#4CAF50"   # green
COL_F1   = "#FF9800"   # orange
DOT_COL  = "#263238"   # dark charcoal for dots
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
# PDF 1: All mutations together — grouped bars + dots
# ============================================================
muts = sorted(results.keys(), key=lambda k: results[k]["f1"], reverse=True)
sens_vals = [results[m]["sensitivity"] for m in muts]
spec_vals = [results[m]["specificity"] for m in muts]
f1_vals   = [results[m]["f1"] for m in muts]

# Summary bars (mean over mutations)
mean_sens = np.mean(sens_vals)
mean_spec = np.mean(spec_vals)
mean_f1   = np.mean(f1_vals)

fig, ax = plt.subplots(figsize=(10, 7))

x = np.arange(3)
bar_width = 0.55
bars = ax.bar(x, [mean_sens, mean_spec, mean_f1], bar_width,
              color=[COL_SENS, COL_SPEC, COL_F1], edgecolor="white", linewidth=1.2,
              alpha=0.85, zorder=3)

# Overlay dots for each mutation
np.random.seed(42)
for i, (vals, col) in enumerate([
    (sens_vals, COL_SENS), (spec_vals, COL_SPEC), (f1_vals, COL_F1)
]):
    jitter = np.random.uniform(-0.15, 0.15, len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals,
               s=28, color=DOT_COL, alpha=0.55, zorder=5, edgecolors="white", linewidths=0.4)

# Labels on bars
for bar, val in zip(bars, [mean_sens, mean_spec, mean_f1]):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
            f"{val:.3f}", ha="center", va="bottom", fontweight="bold", fontsize=12, color="#333")

ax.set_xticks(x)
ax.set_xticklabels(["Sensitivity", "Specificity", "F1 Score"], fontsize=13, fontweight="bold")
ax.set_ylim(0, 1.12)
ax.set_ylabel("Score", fontsize=12)
ax.set_title("BEAT-AML RNA-Seq (No Imputation)\nAll Mutations — Mean ± Per-Mutation Dots",
             fontsize=14, fontweight="bold", pad=15)
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=DOT_COL, markersize=6,
           label=f"Individual mutations (n={len(muts)})", alpha=0.6),
    Line2D([0], [0], color=COL_SENS, lw=8, alpha=0.85, label="Mean Sensitivity"),
    Line2D([0], [0], color=COL_SPEC, lw=8, alpha=0.85, label="Mean Specificity"),
    Line2D([0], [0], color=COL_F1, lw=8, alpha=0.85, label="Mean F1"),
]
ax.legend(handles=legend_elements, loc="lower right", fontsize=9, framealpha=0.9)

# Annotation: n mutations
ax.text(0.01, 0.99, f"BEAT-AML2 cohort (n=707)\n{len(muts)} mutation categories\n5-fold CV-OOF, logL2",
        transform=ax.transAxes, fontsize=8, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#ccc", alpha=0.9))

plt.tight_layout()
pdf1 = os.path.join(OUT_DIR, "beataml_rna_all_mutations.pdf")
with PdfPages(pdf1) as pp:
    pp.savefig(fig, dpi=300)
plt.close(fig)
print(f"Saved: {pdf1}")


# ============================================================
# PDF 2: Per-mutation bar charts (one page per mutation)
# ============================================================
pdf2 = os.path.join(OUT_DIR, "beataml_rna_per_mutation.pdf")
with PdfPages(pdf2) as pp:
    for mut in muts:
        r = results[mut]
        fig, ax = plt.subplots(figsize=(7, 5))

        vals = [r["sensitivity"], r["specificity"], r["f1"]]
        colors = [COL_SENS, COL_SPEC, COL_F1]
        labels = ["Sensitivity", "Specificity", "F1 Score"]
        x = np.arange(3)

        bars = ax.bar(x, vals, 0.55, color=colors, edgecolor="white", linewidth=1.2,
                       alpha=0.88, zorder=3)

        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                    f"{val:.3f}", ha="center", va="bottom", fontweight="bold", fontsize=13,
                    color="#333")

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=12, fontweight="bold")
        ax.set_ylim(0, 1.15)
        ax.set_ylabel("Score", fontsize=11)

        # Clean mutation name for title
        title_name = mut.replace("_", " ").replace("/", " / ")
        ax.set_title(f"BEAT-AML RNA-Seq (No Imputation)\n{title_name}",
                     fontsize=13, fontweight="bold", pad=12)
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))

        # Info box
        info = (f"n+ = {r['n_pos']}   n− = {r['n_neg']}   threshold = {r['threshold']:.3f}\n"
                f"CV AUROC = {r['cv_auroc']}   TP={r['tp']}  FP={r['fp']}  FN={r['fn']}  TN={r['tn']}")
        ax.text(0.5, -0.18, info, transform=ax.transAxes, fontsize=8, ha="center",
                va="top", color="#555",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#ddd", alpha=0.9))

        plt.tight_layout(rect=[0, 0.08, 1, 1])
        pp.savefig(fig, dpi=300)
        plt.close(fig)

print(f"Saved: {pdf2}")
print(f"\nDONE. Charts in: {OUT_DIR}")
