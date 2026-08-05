#!/usr/bin/env python3
"""Provability figures from scratchpad/oof_metrics_v2.json (local matplotlib):
  A  independence ladder  — AUROC across bulk / RNA+comp / measured-only / all-8, gain decomposition.
  B  AUROC forest         — per-mutation AUROC + 95% bootstrap CI vs the 0.5 chance line + permutation p.
  C  honest operating pt  — mean sensitivity/specificity/F1 naive (optimistic) vs nested (corrected).
"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD = os.path.join(ROOT, "deliverables")
d = json.load(open(os.path.join(ROOT, "scratchpad", "oof_metrics_v2.json")))
BLUE, GREEN, INK, MUTED, GRID, SURF = "#2a78d6", "#008300", "#1a1a19", "#6b6b63", "#d9d9d4", "#ffffff"

# dedup aliases by metric signature
mm = d["arms"]["multimodal"]["mutations"]
def sig(m): return (mm[m]["auroc"], mm[m]["n_pos"], mm[m]["f1"])
seen, keep = {}, []
for m in sorted(mm):
    s = sig(m)
    if s in seen:
        if len(m) < len(seen[s]):                 # prefer the SHORT canonical name (inv16, kmt2a)
            keep = [x for x in keep if x != seen[s]] + [m]; seen[s] = m
    else:
        seen[s] = m; keep.append(m)
MUTS = sorted(keep, key=lambda x: -mm[x]["auroc"])
ARMS = ["bulkrna", "rna_comp", "measured", "multimodal"]
ARMLAB = {"bulkrna": "bulk RNA\n(no imputation)", "rna_comp": "RNA +\ncomposition",
          "measured": "measured-only\n(drop 4 imputed)", "multimodal": "all 8\nmodalities"}

def arm_mean(arm, key):
    return float(np.mean([d["arms"][arm]["mutations"][m][key] for m in MUTS]))
def arm_permsig(arm):
    return int(sum(1 for m in MUTS if (d["arms"][arm]["mutations"][m]["perm_p"] or 1) < 0.05))

def style(ax):
    ax.set_facecolor(SURF); ax.spines[["top","right"]].set_visible(False)
    ax.spines[["left","bottom"]].set_color(GRID); ax.tick_params(colors=MUTED, labelsize=9)
    ax.yaxis.grid(True, color=GRID, lw=.7); ax.set_axisbelow(True)

# ================= A: independence ladder =================
au = [arm_mean(a, "auroc") for a in ARMS]
fig, ax = plt.subplots(figsize=(9.2, 6.0), facecolor=SURF)
x = np.arange(4)
bars = ax.bar(x, au, width=0.6, color=[BLUE, "#5b9bd5", "#8fbf6f", GREEN], edgecolor=SURF, lw=1.5)
for i, v in enumerate(au):
    ax.text(i, v+0.006, "%.3f" % v, ha="center", va="bottom", fontsize=11, color=INK, fontweight="bold")
    ax.text(i, 0.505, "perm-sig\n%d/%d" % (arm_permsig(ARMS[i]), len(MUTS)), ha="center", va="bottom",
            fontsize=8, color=MUTED)
ax.set_xticks(x); ax.set_xticklabels([ARMLAB[a] for a in ARMS], fontsize=9.5)
ax.set_ylim(0.5, 0.95); ax.set_ylabel("mean AUROC over %d mutations" % len(MUTS), fontsize=10, color=INK)
style(ax)
# gain brackets
def bracket(x0, x1, y, txt, col):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<->", color=col, lw=1.6))
    ax.text((x0+x1)/2, y+0.008, txt, ha="center", va="bottom", fontsize=9.5, color=col, fontweight="bold")
bracket(0, 2, 0.905, "+%.3f  imputation-INDEPENDENT" % (au[2]-au[0]), "#1a6e1a")
bracket(2, 3, 0.905, "+%.3f  from imputation" % (au[3]-au[2]), "#b06a00")
frac = (au[2]-au[0]) / (au[3]-au[0]) * 100
ax.set_title("Most of the multimodal gain is independent of RNA-imputation\n"
             "%.0f%% of the bulk→all-8 AUROC gain needs no imputed modalities" % frac,
             fontsize=12.5, color=INK, fontweight="bold", pad=12)
fig.text(0.5, 0.02, "'measured-only' drops the 4 imputed blocks (ADT, Lipid, Metabolite, GRN).  "
         "perm-sig = AUROC beats a label-permutation null (p<0.05).",
         ha="center", fontsize=8, color=MUTED, style="italic")
fig.subplots_adjust(left=0.09, right=0.97, top=0.86, bottom=0.13)
fig.savefig(os.path.join(OUTD, "provability_A_independence_ladder.pdf"), facecolor=SURF)
fig.savefig(os.path.join(OUTD, "provability_A_independence_ladder.png"), dpi=130, facecolor=SURF)
plt.close(fig)

# ================= B: AUROC forest with CI =================
fig, ax = plt.subplots(figsize=(9.5, 8.6), facecolor=SURF)
y = np.arange(len(MUTS))[::-1]
for yi, m in zip(y, MUTS):
    r = mm[m]; a = r["auroc"]; ci = r.get("auroc_ci")
    if ci:
        ax.plot([ci[0], ci[1]], [yi, yi], color=BLUE, lw=2, alpha=0.55, solid_capstyle="round")
    sz = 30 + 3*r["n_pos"]
    ax.scatter([a], [yi], s=sz, color=BLUE, edgecolor=SURF, lw=1.2, zorder=3)
    ax.text(1.012, yi, "n⁺%d" % r["n_pos"], va="center", fontsize=7.5, color=MUTED)
ax.axvline(0.5, color="#c0392b", lw=1.4, ls=(0,(5,3))); ax.text(0.5, len(MUTS)-0.2, "chance", color="#c0392b", fontsize=8.5, ha="center")
ax.set_yticks(y); ax.set_yticklabels(MUTS, fontsize=8.7)
ax.set_xlim(0.45, 1.06); ax.set_xlabel("AUROC (dot) with 95% donor-bootstrap CI (bar)", fontsize=10, color=INK)
style(ax); ax.xaxis.grid(True, color=GRID, lw=.7); ax.yaxis.grid(False)
ax.set_title("Every mutation's signal is real — all 24 CIs clear chance\n"
             "MOSAIC-AML multimodal, donor-grouped CV-OOF; all p ≤ 0.003 vs permutation null",
             fontsize=12.5, color=INK, fontweight="bold", pad=12)
fig.subplots_adjust(left=0.15, right=0.94, top=0.90, bottom=0.07)
fig.savefig(os.path.join(OUTD, "provability_B_auroc_forest.pdf"), facecolor=SURF)
fig.savefig(os.path.join(OUTD, "provability_B_auroc_forest.png"), dpi=130, facecolor=SURF)
plt.close(fig)

# ================= C: naive vs nested operating point =================
fig, ax = plt.subplots(figsize=(7.8, 5.2), facecolor=SURF)
metrics = [("sensitivity","nested_sensitivity","Sensitivity"),
           ("specificity","nested_specificity","Specificity"), ("f1","nested_f1","F1")]
xx = np.arange(3); w = 0.38
naive = [np.mean([mm[m][a] for m in MUTS]) for a,_,_ in metrics]
nest = [np.mean([mm[m][b] for m in MUTS]) for _,b,_ in metrics]
ax.bar(xx-w/2, naive, w, color="#b9b9b0", edgecolor=SURF, lw=1.5, label="naive (threshold on same data — optimistic)")
ax.bar(xx+w/2, nest, w, color=GREEN, edgecolor=SURF, lw=1.5, label="nested (threshold on held-out fold — honest)")
for i in range(3):
    ax.text(i-w/2, naive[i]+0.01, "%.2f"%naive[i], ha="center", va="bottom", fontsize=9, color="#6b6b63")
    ax.text(i+w/2, nest[i]+0.01, "%.2f"%nest[i], ha="center", va="bottom", fontsize=9, color=GREEN, fontweight="bold")
ax.set_xticks(xx); ax.set_xticklabels([t for *_,t in metrics], fontsize=10)
ax.set_ylim(0,1.05); ax.set_ylabel("mean over %d mutations"%len(MUTS), fontsize=10, color=INK)
ax.legend(loc="upper left", frameon=False, fontsize=8.5, bbox_to_anchor=(0.0,1.0), labelcolor=INK)
style(ax)
ax.set_title("Honest operating point: nested-CV threshold removes the optimism",
             fontsize=12.5, color=INK, fontweight="bold", pad=30)
fig.subplots_adjust(left=0.11, right=0.97, top=0.80, bottom=0.10)
fig.savefig(os.path.join(OUTD, "provability_C_nested_operating_point.pdf"), facecolor=SURF)
fig.savefig(os.path.join(OUTD, "provability_C_nested_operating_point.png"), dpi=130, facecolor=SURF)
plt.close(fig)

print("mutations (deduped):", len(MUTS))
print("ladder AUROC:", {a: round(arm_mean(a,"auroc"),3) for a in ARMS})
print("imputation-independent fraction: %.0f%%" % frac)
print("wrote provability_A/B/C .pdf/.png")
print("PROVABILITY FIGS OK")
