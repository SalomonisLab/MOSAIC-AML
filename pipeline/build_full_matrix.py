#!/usr/bin/env python3
"""Two models x four cohorts: sensitivity & specificity, overall and per mutation.

  Model 1  Bulk RNA-only predictor  (BeatAML2-trained, no imputation)  -> deliverables/bulk_matrix.json
  Model 2  MOSAIC-AML multimodal    (single-cell atlas, 8 modalities)  -> scratchpad/oof_metrics.json (all-sc)
                                                                          + runs/predict_* (held-out sc)

Cohorts: BeatAML (CV), Leucegene, held-out scRNA (n=29), all scRNA (CV/OOF).
The multimodal model CANNOT run on bulk-only cohorts (BeatAML, Leucegene) — no single-cell modalities.

Writes deliverables/full_matrix.tsv and deliverables/MOSAIC-AML_model_x_cohort_matrix.pdf/.png .
Runs locally (matplotlib 3.10)."""
import os, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)

def canon(c):
    cl = str(c).lower()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl: return "inv16"
    if "kmt2a" in cl: return "KMT2A"
    return str(c).split("_")[0].split("-")[0].upper()

# ---- multimodal all-sc (CV-OOF) ----
om = json.load(open(os.path.join(ROOT, "scratchpad", "oof_metrics.json")))
A = om["arms"]; sh0 = om["paired_summary"]["shared"]
def sig(m): return tuple(A[a]["mutations"][m].get(k) for a in ("bulkrna","multimodal") for k in ("sensitivity","specificity","f1","auroc","n_pos"))
seen, drop = {}, []
for m in sorted(sh0):
    s = sig(m)
    if s in seen: drop.append(m if len(m) <= len(seen[s]) else seen[s]); seen[s] = seen[s] if len(m) <= len(seen[s]) else m
    else: seen[s] = m
shared = [m for m in sh0 if m not in set(drop)]
mm_all = {canon(m): A["multimodal"]["mutations"][m] for m in shared}

# ---- multimodal held-out sc (from reports) ----
conf = {}
for f in sorted(glob.glob(os.path.join(ROOT, "runs", "predict_*", "patient_report.json"))):
    r = json.load(open(f))
    for p in (r.get("mutation_predictions") or []):
        tl = p.get("true_label")
        if tl is None: continue
        d = conf.setdefault(canon(p["mutation"]), [0, 0, 0, 0])
        pos = (tl == "present"); cp = (p.get("call") == "present")
        d[0 if (pos and cp) else 2 if (pos and not cp) else 1 if cp else 3] += 1
mm_ho = {}
for m, (tp, fp, fn, tn) in conf.items():
    if tp + fn >= 1:
        mm_ho[m] = {"sensitivity": tp/(tp+fn), "specificity": (tn/(tn+fp) if tn+fp else float("nan")),
                    "n_pos": tp+fn, "tp": tp, "fp": fp, "fn": fn, "tn": tn}

# ---- bulk model, 4 cohorts ----
bm = json.load(open(os.path.join(OUTD, "bulk_matrix.json")))["cohorts"]
def bulk_percat(coh):
    return {canon(c): v for c, v in bm[coh]["per_category"].items()}
CELLS = [("Bulk·BeatAML(CV)", bulk_percat("BeatAML_CV"), bm["BeatAML_CV"]["overall"]),
         ("Bulk·Leucegene", bulk_percat("Leucegene"), bm["Leucegene"]["overall"]),
         ("Bulk·held-out sc", bulk_percat("heldout_scRNA"), bm["heldout_scRNA"]["overall"]),
         ("Bulk·all sc", bulk_percat("all_scRNA"), bm["all_scRNA"]["overall"]),
         ("MM·held-out sc", mm_ho, None),
         ("MM·all sc(CV)", mm_all, None)]

def overall(cellmap, stored):
    if stored: return stored["mean_sensitivity"], stored["mean_specificity"], stored.get("mean_auroc")
    se = np.nanmean([v["sensitivity"] for v in cellmap.values()])
    sp = np.nanmean([v["specificity"] for v in cellmap.values()])
    au = [v.get("auroc") for v in cellmap.values() if v.get("auroc") is not None]
    return se, sp, (np.mean(au) if au else None)

# ---------------- combined TSV ----------------
KEY = ["inv16","NPM1","FLT3-ITD","FLT3","TP53","del7","del5","trisomy8","complex","KMT2A",
       "TET2","DNMT3A","ASXL1","NRAS","KRAS","RUNX1","IDH1","IDH2","NPM1","WT1","CEBPA","KIT","GATA2","PTPN11"]
KEY = list(dict.fromkeys(KEY))
with open(os.path.join(OUTD, "full_matrix.tsv"), "w") as fh:
    fh.write("cell\tscope\tmutation\tsensitivity\tspecificity\tn_pos\n")
    for name, cm, st in CELLS:
        se, sp, au = overall(cm, st)
        fh.write("%s\tOVERALL\t(mean)\t%.3f\t%.3f\t\n" % (name, se, sp))
        for mut in KEY:
            v = cm.get(mut)
            if v:
                fh.write("%s\tper-mut\t%s\t%.3f\t%.3f\t%s\n" %
                         (name, mut, v["sensitivity"], v["specificity"], v.get("n_pos","")))

# ---------------- figure ----------------
C_SENS, C_SPEC = "#2a78d6", "#008300"         # validated categorical slots 1 & 2
INK, MUTED, GRID, SURF, NA = "#1a1a19", "#6b6b63", "#d9d9d4", "#ffffff", "#eceae2"
fig = plt.figure(figsize=(13.5, 8.6), facecolor=SURF)
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.25], hspace=0.42)

# -- panel A: overall sens/spec per cell --
axA = fig.add_subplot(gs[0]); axA.set_facecolor(SURF)
names = [c[0] for c in CELLS]; x = np.arange(len(names)); w = 0.38
se = [overall(c[1], c[2])[0] for c in CELLS]; sp = [overall(c[1], c[2])[1] for c in CELLS]
au = [overall(c[1], c[2])[2] for c in CELLS]
axA.bar(x - w/2, se, w, color=C_SENS, edgecolor=SURF, lw=1.5, label="sensitivity")
axA.bar(x + w/2, sp, w, color=C_SPEC, edgecolor=SURF, lw=1.5, label="specificity")
for i in range(len(names)):
    axA.text(x[i]-w/2, se[i]+0.02, "%.2f"%se[i], ha="center", va="bottom", fontsize=8, color=C_SENS, fontweight="bold")
    axA.text(x[i]+w/2, sp[i]+0.02, "%.2f"%sp[i], ha="center", va="bottom", fontsize=8, color=C_SPEC, fontweight="bold")
    if au[i] is not None:
        axA.text(x[i], -0.13, "AUROC\n%.2f"%au[i], ha="center", va="top", fontsize=7.5, color=MUTED)
# N/A columns for multimodal on bulk cohorts
for xi, lab in [(len(names)+0.0,"MM·BeatAML"),(len(names)+1.0,"MM·Leucegene")]:
    axA.axvspan(xi-0.5, xi+0.5, color=NA, zorder=0)
    axA.text(xi, 0.5, "N/A\nbulk cohort:\nno single-cell\nmodalities", ha="center", va="center",
             fontsize=8, color=MUTED, style="italic")
axA.set_xlim(-0.6, len(names)+1.6); axA.set_ylim(0, 1.06)
axA.set_xticks(list(x)+[len(names),len(names)+1]); axA.set_xticklabels(names+["MM·BeatAML","MM·Leucegene"], fontsize=8.7)
axA.set_ylabel("score", fontsize=10, color=INK); axA.set_yticks([0,.25,.5,.75,1.0])
axA.spines[["top","right"]].set_visible(False); axA.spines[["left","bottom"]].set_color(GRID)
axA.tick_params(colors=MUTED, labelsize=8.5); axA.yaxis.grid(True, color=GRID, lw=.7); axA.set_axisbelow(True)
axA.axvline(3.5, color=GRID, lw=1, ls=(0,(4,3)))
axA.text(1.5, 1.10, "Bulk RNA-only model ①  (BeatAML-trained)", ha="center", fontsize=9.5, color=INK, fontweight="bold", transform=axA.transData)
axA.text(5.4, 1.10, "MOSAIC-AML multimodal ②  (single-cell)", ha="center", fontsize=9.5, color=INK, fontweight="bold", transform=axA.transData)
handlesA = [Patch(facecolor=C_SENS, label="sensitivity"), Patch(facecolor=C_SPEC, label="specificity")]
fig.legend(handles=handlesA, loc="upper center", ncol=2, frameon=False, fontsize=9.5,
           bbox_to_anchor=(0.5, 0.945), labelcolor=INK)

# -- panel B: per-mutation sensitivity heatmap --
axB = fig.add_subplot(gs[1]); axB.set_facecolor(SURF)
muts = [m for m in KEY if any((c[1].get(m)) for c in CELLS)]
M = np.full((len(muts), len(CELLS)), np.nan)
for j, (nm, cm, st) in enumerate(CELLS):
    for i, mut in enumerate(muts):
        v = cm.get(mut)
        if v: M[i, j] = v["sensitivity"]
im = axB.imshow(M, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
axB.set_xticks(range(len(CELLS))); axB.set_xticklabels(names, fontsize=8.5, rotation=25, ha="right")
axB.set_yticks(range(len(muts))); axB.set_yticklabels(muts, fontsize=8.5)
for i in range(len(muts)):
    for j in range(len(CELLS)):
        if not np.isnan(M[i, j]):
            axB.text(j, i, "%.2f"%M[i,j], ha="center", va="center", fontsize=7,
                     color="white" if M[i,j] > 0.55 else INK)
        else:
            axB.add_patch(plt.Rectangle((j-0.5,i-0.5),1,1,color=NA,zorder=1))
axB.set_title("Per-mutation SENSITIVITY  (blank = too few positives to evaluate in that cohort)",
              fontsize=11.5, color=INK, fontweight="bold", pad=8)
cb = fig.colorbar(im, ax=axB, fraction=0.025, pad=0.01); cb.set_label("sensitivity", fontsize=8, color=MUTED)
cb.ax.tick_params(labelsize=7, colors=MUTED)

fig.suptitle("MOSAIC-AML — mutation-calling performance by model and cohort", fontsize=15.5,
             color=INK, fontweight="bold", y=0.988)
fig.text(0.5, 0.008,
         "Bulk model ① is BeatAML2-trained (bulk RNA, no imputation); multimodal ② is single-cell-trained (8 modalities). "
         "Both operate at the F1-max threshold. The multimodal model has no bulk-cohort columns — bulk cohorts lack single-cell modalities.",
         ha="center", fontsize=7.6, color=MUTED, style="italic")
fig.subplots_adjust(left=0.10, right=0.98, top=0.86, bottom=0.11)
for ext, dpi in [("pdf", None), ("png", 130)]:
    fig.savefig(os.path.join(OUTD, "MOSAIC-AML_model_x_cohort_matrix.%s" % ext),
                **({"dpi": dpi} if dpi else {}), facecolor=SURF)
print("wrote deliverables/MOSAIC-AML_model_x_cohort_matrix.pdf/.png + full_matrix.tsv")
print("MATRIX OK")
