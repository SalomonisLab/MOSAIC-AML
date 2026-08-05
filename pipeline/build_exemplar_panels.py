#!/usr/bin/env python3
"""One chart per exemplar driver: how every modality performs for that lesion.

Each panel ranks all nine measurement blocks by their STANDALONE donor-grouped CV-OOF AUROC for that
driver, marks the winner, and prints the weight the deployed fusion assigns each block — so the two
independent lines of evidence (does it predict alone? does the model rely on it?) are visible together.

  python build_exemplar_panels.py -> deliverables/figures/Fig7_exemplar_modality_profiles.pdf/.png
"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "deliverables"); FIGD = os.path.join(D, "figures")
MM = 1/25.4
INK, MUTED, GRID = "#1a1a19", "#6b6b63", "#d9d9d4"
WIN, OTH, WT = "#c0392b", "#9ec9e2", "#1a6e1a"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
    "font.size":7,"axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
    "legend.fontsize":6.5,"axes.linewidth":0.6,"pdf.fonttype":42,"savefig.dpi":600})

B = json.load(open(os.path.join(D, "modality_breakdown_current.json")))["drivers"]
P = json.load(open(os.path.join(D, "production_fused_model.json")))["per_mutation"]
LAB = {"RNA":"RNA (cell-state)","BulkRNA":"RNA (whole-sample)","Composition":"cell-state composition",
       "ADT":"ADT (surface protein)","Lipid":"lipid","Metabolite":"metabolite","GRN":"GRN",
       "LSC":"LSC score","Cell-comm":"cell–cell communication"}

# exemplars chosen in build_modality_figure.py: best driver per modality by margin over runner-up
EX = [("GATA2","ADT"), ("FLT3-TKD","RNA"), ("RUNX1","Metabolite"),
      ("TET2","GRN"), ("ASXL1","Cell-comm"), ("SRSF2","Lipid")]

def fused(m): return ((P.get(m,{}) or {}).get("fused_all") or {}).get("auroc")
def npos(m):  return ((P.get(m,{}) or {}).get("fused_all") or {}).get("n_pos_atlas")

fig, axes = plt.subplots(2, 3, figsize=(180*MM, 122*MM))
axes = axes.ravel()
for ax, (drv, mod) in zip(axes, EX):
    r = B[drv]; sa = r["standalone_auroc"]; w = r.get("fusion_weights") or {}
    items = sorted(sa.items(), key=lambda kv: kv[1])          # ascending -> best at top
    y = np.arange(len(items))
    vals = [v for _, v in items]
    cols = [WIN if k == mod else OTH for k, _ in items]
    ax.barh(y, [v-0.5 for v in vals], left=0.5, height=0.68, color=cols, edgecolor="white", lw=0.4)
    for i, (k, v) in enumerate(items):
        ax.text(v+0.006, i, "%.2f" % v, va="center", fontsize=5.6,
                color=WIN if k == mod else MUTED, fontweight="bold" if k == mod else "normal")
        if w.get(k, 0) > 0.01:                                 # fusion weight, where the model uses it
            ax.text(0.508, i, "w %.2f" % w[k], va="center", fontsize=5.0, color=WT)
    ax.set_yticks(y); ax.set_yticklabels([LAB[k] for k, _ in items], fontsize=5.8)
    ax.set_xlim(0.5, 1.06); ax.set_xticks([0.5,0.6,0.7,0.8,0.9,1.0])
    ax.axvline(0.5, color=GRID, lw=0.6)
    best = items[-1][1]; second = items[-2][1]
    ax.set_title("%s   (n+ %d,  fused %.3f)\nbest: %s  %.3f  (+%.3f)"
                 % (drv, npos(drv) or 0, fused(drv) or 0, LAB[mod], best, best-second),
                 fontsize=6.4, color=INK, fontweight="bold", loc="left", pad=4)
    ax.set_facecolor("white")
    for s in ("top","right"): ax.spines[s].set_visible(False)
    for s in ("left","bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED); ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
for ax in axes[len(EX):]: ax.axis("off")
fig.supxlabel("standalone AUROC of each modality alone (donor-grouped CV-OOF)", fontsize=7.5, y=0.028)
fig.legend(handles=[Line2D([],[],marker="s",ls="",color=WIN,label="best modality for this driver",ms=5),
                    Line2D([],[],marker="s",ls="",color=OTH,label="other modalities",ms=5),
                    Line2D([],[],marker="",ls="",color=WT,label="w = weight in the deployed fusion")],
           frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.005))
fig.subplots_adjust(left=0.115, right=0.985, top=0.93, bottom=0.105, hspace=0.70, wspace=0.75)
for ext, dpi in (("pdf", None), ("png", 200)):
    fig.savefig(os.path.join(FIGD, "Fig7_exemplar_modality_profiles.%s" % ext),
                **({"dpi": dpi} if dpi else {}), bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print("wrote Fig7_exemplar_modality_profiles.pdf/.png")
for drv, mod in EX:
    sa = B[drv]["standalone_auroc"]; s = sorted(sa.items(), key=lambda kv: -kv[1])
    print("  %-10s winner %-12s %.3f  (runner-up %s %.3f)  fused %.3f" %
          (drv, mod, sa[mod], s[1][0], s[1][1], fused(drv) or 0))
