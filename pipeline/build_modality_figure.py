#!/usr/bin/env python3
"""Figure + table: which measurement type carries which driver.

a  heatmap of standalone donor-grouped CV-OOF AUROC for every (driver x modality)
b  exemplar drivers, one per modality, showing the full modality profile and the winning block

Also writes the complete matrix as a TSV for supplementary submission.
"""
import os, json, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "deliverables")
FIGD = os.path.join(D, "figures"); TABD = os.path.join(D, "tables")
os.makedirs(FIGD, exist_ok=True); os.makedirs(TABD, exist_ok=True)
MM = 1/25.4
INK, MUTED, GRID = "#1a1a19", "#6b6b63", "#d9d9d4"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
    "font.size":7,"axes.labelsize":7.5,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
    "legend.fontsize":6.5,"axes.linewidth":0.6,"pdf.fonttype":42,"savefig.dpi":600})

B = json.load(open(os.path.join(D, "modality_breakdown_current.json")))
P = json.load(open(os.path.join(D, "production_fused_model.json")))["per_mutation"]
MODS = ["RNA","BulkRNA","Composition","ADT","Lipid","Metabolite","GRN","LSC","Cell-comm"]
LAB  = {"RNA":"RNA\n(cell-state)","BulkRNA":"RNA\n(whole-sample)","Composition":"cell-state\ncomposition",
        "ADT":"ADT\n(surface)","Lipid":"lipid","Metabolite":"metabolite","GRN":"GRN","LSC":"LSC",
        "Cell-comm":"cell–cell\ncomm."}

drivers = B["drivers"]
def fused(m): return ((P.get(m,{}) or {}).get("fused_all") or {}).get("auroc")
# de-duplicate alias lesions, keep the shorter canonical name
seen, keep = {}, []
for m, r in drivers.items():
    k = (r["best_single_auroc"], r["n_pos"])
    if k in seen:
        if len(m) < len(seen[k]):
            keep = [x for x in keep if x != seen[k]] + [m]; seen[k] = m
    else:
        seen[k] = m; keep.append(m)
order = sorted(keep, key=lambda m: -(fused(m) or 0))

# ------------------------------------------------ TSV (full matrix)
with open(os.path.join(TABD, "SuppTable5_modality_breakdown.tsv"), "w", encoding="utf-8") as fh:
    fh.write("driver\tn_positive\tfused_AUROC\tbest_modality\tbest_single_AUROC\tmargin_over_runner_up\t"
             + "\t".join(MODS) + "\t" + "\t".join("weight_"+m for m in MODS) + "\n")
    for m in order:
        r = drivers[m]; sa = r["standalone_auroc"]; w = r.get("fusion_weights") or {}
        s = sorted(sa.items(), key=lambda kv: -kv[1])
        margin = (s[0][1]-s[1][1]) if len(s) > 1 else 0
        fh.write("\t".join([m, str(r["n_pos"]), "%.4f" % (fused(m) or 0), r["best_modality"],
                            "%.4f" % r["best_single_auroc"], "%+.4f" % margin]
                 + ["%.4f" % sa[x] if x in sa else "" for x in MODS]
                 + ["%.3f" % w[x] if x in w else "" for x in MODS]) + "\n")
print("  wrote SuppTable5_modality_breakdown.tsv (%d drivers)" % len(order))

# ------------------------------------------------ exemplars: best driver per modality by margin
best_for = {}
for m in order:
    sa = drivers[m]["standalone_auroc"]
    s = sorted(sa.items(), key=lambda kv: -kv[1])
    if len(s) < 2: continue
    mod, a, margin = s[0][0], s[0][1], s[0][1]-s[1][1]
    if mod not in best_for or margin > best_for[mod][2]:
        best_for[mod] = (m, a, margin, s[1][0], s[1][1])
EX_ORDER = ["ADT","RNA","Lipid","Metabolite","Cell-comm","GRN"]
EX = [(mod,)+best_for[mod] for mod in EX_ORDER if mod in best_for]

# ------------------------------------------------ figure
fig = plt.figure(figsize=(180*MM, 150*MM))
gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.0], hspace=0.30)

# a: heatmap
axA = fig.add_subplot(gs[0])
Mx = np.full((len(order), len(MODS)), np.nan)
for i, m in enumerate(order):
    for j, mod in enumerate(MODS):
        v = drivers[m]["standalone_auroc"].get(mod)
        if v is not None: Mx[i, j] = v
im = axA.imshow(Mx, aspect="auto", cmap="YlGnBu", vmin=0.5, vmax=1.0)
axA.set_xticks(range(len(MODS))); axA.set_xticklabels([LAB[m] for m in MODS], fontsize=6)
axA.set_yticks(range(len(order))); axA.set_yticklabels(order, fontsize=5.6)
for i, m in enumerate(order):                                   # ring the winning modality
    bm = drivers[m]["best_modality"]
    if bm in MODS:
        j = MODS.index(bm)
        axA.add_patch(plt.Rectangle((j-.5, i-.5), 1, 1, fill=False, ec="#c0392b", lw=1.1))
for i in range(len(order)):
    for j in range(len(MODS)):
        if not np.isnan(Mx[i, j]):
            axA.text(j, i, "%.2f" % Mx[i, j], ha="center", va="center", fontsize=4.2,
                     color="white" if Mx[i, j] > 0.85 else INK)
cb = fig.colorbar(im, ax=axA, fraction=0.018, pad=0.01)
cb.set_label("standalone AUROC", fontsize=6.5); cb.ax.tick_params(labelsize=5.5)
axA.set_title("a   Standalone performance of each modality for each driver "
              "(red outline = best single modality)", fontsize=8, fontweight="bold",
              color=INK, loc="left", pad=6)

# b: exemplars
axB = fig.add_subplot(gs[1])
n = len(EX); w = 0.82/len(MODS)
for i, (mod, drv, a, margin, r2, a2) in enumerate(EX):
    sa = drivers[drv]["standalone_auroc"]
    for j, mm2 in enumerate(MODS):
        v = sa.get(mm2)
        if v is None: continue
        col = "#c0392b" if mm2 == mod else "#9ec9e2"
        xpos = i - 0.41 + j*w + w/2
        axB.bar(xpos, v-0.5, width=w*0.9, bottom=0.5, color=col, edgecolor="white", lw=0.3)
        if mm2 == mod:                                   # annotate only the winner, above its own bar
            axB.text(xpos, v+0.008, "%.2f" % v, ha="center", va="bottom",
                     fontsize=5.6, color="#c0392b", fontweight="bold")
    if i:                                                # separator between exemplars
        axB.axvline(i-0.5, color=GRID, lw=0.5)
axB.set_xticks(range(n))
axB.set_xticklabels(["%s\n%s  (+%.3f)" % (e[1], LAB[e[0]].replace("\n", " "), e[3]) for e in EX],
                    fontsize=6.2)
axB.set_ylim(0.5, 1.0); axB.set_ylabel("standalone AUROC", fontsize=7.5)
axB.set_xlim(-0.55, n-0.45)
axB.set_facecolor("white")
for s in ("top","right"): axB.spines[s].set_visible(False)
for s in ("left","bottom"): axB.spines[s].set_color(GRID)
axB.tick_params(colors=MUTED); axB.yaxis.grid(True, color=GRID, lw=.45); axB.set_axisbelow(True)
axB.legend(handles=[Line2D([],[],marker="s",ls="",color="#c0392b",label="winning modality",ms=5),
                    Line2D([],[],marker="s",ls="",color="#9ec9e2",label="other modalities",ms=5)],
           frameon=False, loc="lower right", ncol=2, handletextpad=.4)
axB.set_title("b   Exemplar drivers — each carried by a different measurement type "
              "(bars = all nine modalities in the order of a; margin over runner-up in brackets)",
              fontsize=8, fontweight="bold", color=INK, loc="left", pad=6)
fig.subplots_adjust(left=0.135, right=0.965, top=0.955, bottom=0.055)
for ext, dpi in (("pdf", None), ("png", 200)):
    fig.savefig(os.path.join(FIGD, "Fig6_modality_specificity.%s" % ext),
                **({"dpi": dpi} if dpi else {}), bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print("  wrote Fig6_modality_specificity.pdf/.png")
print("\nexemplars:")
for mod, drv, a, margin, r2, a2 in EX:
    print("   %-12s %-10s standalone %.3f  (+%.3f over %s)  fused %.3f" %
          (mod, drv, a, margin, r2, fused(drv) or 0))
never = [m for m in MODS if m not in best_for]
print("\nnever the best single modality:", never or "none")
