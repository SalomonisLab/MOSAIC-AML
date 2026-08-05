#!/usr/bin/env python3
"""Figure for the BeatAML imputation experiment, scored by mean per-mutation ACCURACY (at the F1-max
operating point) rather than AUROC. Shows the trivial 'call-all-absent' baseline that dominates accuracy
on rare mutations. (local matplotlib)"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD = os.path.join(ROOT, "deliverables")
d = json.load(open(os.path.join(OUTD, "beataml_impute_experiment.json")))
m = d["arm_means"]; pm = d["per_mutation"]
prev = [pm[k]["rna_only"]["prevalence"] for k in pm]
BASE = float(np.mean([max(p, 1 - p) for p in prev]))          # 'call the majority class' accuracy, averaged
GREY, GREEN, RED, INK, MUTED, GRID, SURF = "#b9b9b0", "#008300", "#c0392b", "#1a1a19", "#6b6b63", "#d9d9d4", "#ffffff"

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.6), facecolor=SURF, gridspec_kw={"width_ratios": [1.15, 1]})

# ---- A: three-arm accuracy bars + trivial baseline ----
labs = ["RNA only", "RNA + imputed\n(ADT·Lipid·Metab·GRN)", "RNA + random\nnonlinear (control)"]
vals = [m["rna_only"]["accuracy"], m["rna_imputed"]["accuracy"], m["rna_random"]["accuracy"]]
bars = axA.bar([0, 1, 2], vals, width=0.62, color=[GREY, GREEN, GREEN], edgecolor=SURF, lw=1.5)
bars[2].set_hatch("////"); bars[2].set_edgecolor("#ffffff")
for i, v in enumerate(vals):
    axA.text(i, v + 0.0015, "%.3f" % v, ha="center", va="bottom", fontsize=12, color=INK, fontweight="bold")
axA.axhline(BASE, color=RED, lw=1.6, ls=(0, (5, 3)))
axA.text(2.46, BASE, "call-all-absent\nbaseline %.3f" % BASE, color=RED, fontsize=9, va="center", ha="left", fontweight="bold")
axA.set_xticks([0, 1, 2]); axA.set_xticklabels(labs, fontsize=9.5)
axA.set_ylim(0.88, 0.965); axA.set_ylabel("mean per-mutation accuracy (47 mut, F1-max operating point)", fontsize=9.5, color=INK)
axA.set_xlim(-0.6, 3.15)
axA.set_facecolor(SURF); axA.spines[["top", "right"]].set_visible(False)
axA.spines[["left", "bottom"]].set_color(GRID); axA.tick_params(colors=MUTED, labelsize=9)
axA.yaxis.grid(True, color=GRID, lw=.7); axA.set_axisbelow(True)
def bracket(x0, x1, y, txt, col):
    axA.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<->", color=col, lw=1.5))
    axA.text((x0 + x1) / 2, y + 0.0016, txt, ha="center", va="bottom", fontsize=9, color=col, fontweight="bold")
bracket(0, 1, 0.947, "+%.3f  p=%.3f" % (d["imputed_vs_only"]["delta"], d["imputed_vs_only"]["p"]), "#1a6e1a")
axA.text(1.5, 0.958, "imputed ≈ random   Δ%+.3f  p=%.2f (n.s.)" % (d["imputed_vs_random_REAL_BIOLOGY"]["delta"], d["imputed_vs_random_REAL_BIOLOGY"]["p"]),
         ha="center", va="center", fontsize=9, color="#b06a00", fontweight="bold")
axA.set_title("Scored by accuracy: same story, and the do-nothing baseline wins", fontsize=12, color=INK, fontweight="bold", pad=10)

# ---- B: per-mutation imputed vs random accuracy ----
xs = [pm[k]["rna_random"]["accuracy"] for k in pm]
ys = [pm[k]["rna_imputed"]["accuracy"] for k in pm]
lo = min(min(xs), min(ys)) - 0.02
axB.plot([lo, 1.001], [lo, 1.001], color=MUTED, ls=(0, (5, 4)), lw=1.2)
axB.scatter(xs, ys, s=46, color=GREEN, edgecolor=SURF, lw=1, alpha=0.85)
axB.set_xlim(lo, 1.005); axB.set_ylim(lo, 1.005)
axB.set_xlabel("accuracy · RNA + random nonlinear", fontsize=10); axB.set_ylabel("accuracy · RNA + imputed", fontsize=10)
axB.set_facecolor(SURF); axB.spines[["top", "right"]].set_visible(False)
axB.spines[["left", "bottom"]].set_color(GRID); axB.tick_params(colors=MUTED, labelsize=9)
axB.grid(True, color=GRID, lw=.6); axB.set_axisbelow(True)
axB.set_title("Per mutation: imputed vs random accuracy", fontsize=12, color=INK, fontweight="bold", pad=10)

fig.suptitle("BeatAML imputation experiment — mean per-mutation ACCURACY (not AUROC)", fontsize=13.5, color=INK, fontweight="bold", y=0.985)
fig.text(0.5, 0.012, "Accuracy is prevalence-dominated (mean prevalence %.1f%%): the trivial call-all-absent classifier scores %.3f, "
         "beating every model at its operating point — so accuracy barely discriminates here. Still, the ordering matches AUROC: imputed ≈ random "
         "(Δ%+.3f, p=%.2f) → nonlinear feature-expansion, not biology." % (100*np.mean(prev), BASE,
         d["imputed_vs_random_REAL_BIOLOGY"]["delta"], d["imputed_vs_random_REAL_BIOLOGY"]["p"]),
         ha="center", fontsize=8, color=MUTED, style="italic", wrap=True)
fig.subplots_adjust(left=0.075, right=0.98, top=0.85, bottom=0.16, wspace=0.22)
for ext, dpi in [("pdf", None), ("png", 130)]:
    fig.savefig(os.path.join(OUTD, "beataml_impute_experiment.%s" % ext), **({"dpi": dpi} if dpi else {}), facecolor=SURF)
print("accuracy: only %.4f imputed %.4f random %.4f | baseline %.4f" % (vals[0], vals[1], vals[2], BASE))
print("wrote deliverables/beataml_impute_experiment.pdf/.png")
